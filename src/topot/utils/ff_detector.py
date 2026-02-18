"""Force field detection from topology files with multiple fallback strategies"""

from pathlib import Path
import re


def detect_force_field(top_file, ff_dir):
    """
    Detect force field from topology with multiple strategies:
    1. Parse #include directives (primary)
    2. Try similar FF names if exact match fails (fallback 1)
    3. Infer from atom types in topology (fallback 2)
    4. Return 'unknown' with diagnostic info

    ff_dir can be either:
    - A directory containing .ff subdirectories (standard case)
    - A .ff directory itself (if user passes a .ff dir directly)

    Returns dict with 'name', 'path', 'mutres_path', 'detection_method'
    """
    ff_info = {
        'name': 'unknown',
        'path': None,
        'mutres_path': None,
        'detection_method': None,
        'available_ff': []
    }

    ff_dir = Path(ff_dir)

    # If ff_dir is a .ff directory itself, use its parent
    if ff_dir.name.endswith('.ff'):
        ff_search_dir = ff_dir.parent
        ff_dir_name = ff_dir.name
    else:
        ff_search_dir = ff_dir
        ff_dir_name = None

    # Collect available FFs for fallback
    available_ffs = list(ff_search_dir.glob('*.ff'))
    ff_info['available_ff'] = [f.name.replace('.ff', '') for f in available_ffs]

    # Read topology file
    with open(top_file) as f:
        content = f.read()

    # METHOD 1: Parse #include directives
    result = _detect_from_includes(content, ff_search_dir, ff_dir_name)
    if result['found']:
        ff_info.update(result)
        ff_info['detection_method'] = '#include directive'
        return ff_info

    # METHOD 2: Try similar FF names (fallback 1)
    # Example: amber99sb-star-ildn-mut-charged -> try amber99sb-star-ildn-mut
    if result.get('requested_name'):
        result = _detect_from_similar_names(result['requested_name'], ff_search_dir)
        if result['found']:
            ff_info.update(result)
            ff_info['detection_method'] = 'similar name (suffix stripped)'
            return ff_info

    # METHOD 3: Infer from atom types in topology
    result = _detect_from_atom_types(content, ff_search_dir)
    if result['found']:
        ff_info.update(result)
        ff_info['detection_method'] = 'atom type inference'
        return ff_info

    # No FF found
    ff_info['detection_method'] = 'none (returned unknown)'
    return ff_info


def _detect_from_includes(content, ff_search_dir, ff_dir_name=None):
    """Strategy 1: Parse #include directives (primary method)"""
    result = {'found': False, 'requested_name': None, 'name': None, 'path': None, 'mutres_path': None}

    # Look for force field includes
    includes = re.findall(r'#include\s*["\']([^"\']+)["\']', content)

    for include in includes:
        # Extract force field name
        if '.ff/' in include:
            # Extract the part before .ff/
            ff_part = include.split('.ff/')[0]
            # Get just the directory name (last component)
            ff_name = Path(ff_part).name + '.ff' if '/' in ff_part or '\\' in ff_part else ff_part + '.ff'

            result['requested_name'] = ff_name.replace('.ff', '')

            # If we know the .ff dir, check if it matches
            if ff_dir_name and ff_name != ff_dir_name:
                continue

            # Search for this FF in ff_search_dir
            ff_path = ff_search_dir / ff_name
            if ff_path.exists():
                result['found'] = True
                result['name'] = ff_name.replace('.ff', '')
                result['path'] = str(ff_path)

                # Look for mutres.mtp
                mutres_path = ff_path / 'mutres.mtp'
                if mutres_path.exists():
                    result['mutres_path'] = str(mutres_path)
                break

    return result


def _detect_from_similar_names(requested_ff_name, ff_search_dir):
    """Strategy 2: Try similar FF names by stripping suffixes"""
    result = {'found': False, 'name': None, 'path': None, 'mutres_path': None}

    # List of common suffixes to try removing
    suffixes_to_try = [
        '-charged',
        '-mut-charged',
        '-dna',
        '-bsc1',
        '-gaff',
    ]

    # Try progressively removing suffixes
    base_name = requested_ff_name
    for suffix in suffixes_to_try:
        if suffix in base_name:
            base_name = base_name.replace(suffix, '')
            ff_name = base_name + '.ff'
            ff_path = ff_search_dir / ff_name

            if ff_path.exists():
                result['found'] = True
                result['name'] = base_name
                result['path'] = str(ff_path)

                # Look for mutres.mtp
                mutres_path = ff_path / 'mutres.mtp'
                if mutres_path.exists():
                    result['mutres_path'] = str(mutres_path)
                break

    return result


def _detect_from_atom_types(content, ff_search_dir):
    """Strategy 3: Infer FF from atom types in topology"""
    result = {'found': False, 'name': None, 'path': None, 'mutres_path': None}

    # Extract atom types from [ atoms ] section
    atom_types = _extract_atom_types(content)
    if not atom_types:
        return result

    # Check each available FF to see which one is most compatible
    for ff_path in sorted(ff_search_dir.glob('*.ff')):
        ff_atom_types = _get_ff_atom_types(ff_path)
        if ff_atom_types:
            # Check overlap: how many topology atom types exist in this FF?
            overlap = len(atom_types & ff_atom_types)
            if overlap > len(atom_types) * 0.5:  # At least 50% match
                result['found'] = True
                result['name'] = ff_path.name.replace('.ff', '')
                result['path'] = str(ff_path)

                # Look for mutres.mtp
                mutres_path = ff_path / 'mutres.mtp'
                if mutres_path.exists():
                    result['mutres_path'] = str(mutres_path)
                break

    return result


def _extract_atom_types(content):
    """Extract atom types from topology [ atoms ] section"""
    atom_types = set()

    # Find the [ atoms ] section
    atoms_match = re.search(r'\[\s*atoms\s*\](.*?)(?=\[|$)', content, re.DOTALL | re.IGNORECASE)
    if not atoms_match:
        return atom_types

    atoms_section = atoms_match.group(1)

    # Parse atom type lines (column 2)
    # Format: nr type resnr residue atom cgnr charge mass [typeB chargeB massB]
    for line in atoms_section.split('\n'):
        # Skip comments and empty lines
        if line.strip().startswith(';') or not line.strip():
            continue

        parts = line.split()
        if len(parts) >= 2:
            # Type is second column
            atom_type = parts[1]
            # Skip DUM_ atoms as they're not real
            if not atom_type.startswith('DUM_'):
                atom_types.add(atom_type)

    return atom_types


def _get_ff_atom_types(ff_path):
    """Get atom types defined in a force field"""
    atom_types = set()

    # Check for atomtypes.atp file
    atp_file = ff_path / 'atomtypes.atp'
    if atp_file.exists():
        try:
            with open(atp_file) as f:
                for line in f:
                    if line.strip().startswith(';') or not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) >= 1:
                        atom_types.add(parts[0])
        except Exception:
            pass

    return atom_types
