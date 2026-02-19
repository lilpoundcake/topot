"""Parse GROMACS topology files to extract atom information"""

from pathlib import Path
import re


def parse_topology(top_file):
    """
    Parse topology file and extract atoms with dual topology info.

    Each .itp file may restart atom numbering from 1, so we apply an offset
    to ensure unique keys when combining atoms from multiple files.

    Returns dict mapping global_atom_index -> {
        'index': int,
        'type': str,           # State A type
        'typeB': str | None,   # State B type
        'resnr': int,
        'resname': str,
        'atomname': str,
        'charge': float,
        'chargeB': float | None,
        'mass': float,
        'massB': float | None
    }
    """
    atoms = {}

    # Collect all .itp files referenced in the topology
    itp_files = [top_file]
    with open(top_file) as f:
        content = f.read()

    # Find all #include directives
    includes = re.findall(r'#include\s*["\']([^"\']+)["\']', content)
    top_dir = Path(top_file).parent
    for inc in includes:
        inc_path = top_dir / inc
        if inc_path.exists():
            itp_files.append(inc_path)

    # Parse each file to find [ atoms ] sections
    # Use offset to handle .itp files that restart atom numbering from 1
    # Track chain_idx so atoms can be associated with their source chain
    offset = 0
    chain_idx = 0
    for itp_file in itp_files:
        atoms_section = extract_section(itp_file, '[ atoms ]')
        if atoms_section:
            parsed = parse_atoms_section(atoms_section)
            if parsed:
                # Apply offset to atom indices to avoid key collisions
                for orig_idx, atom_data in parsed.items():
                    new_idx = orig_idx + offset
                    atom_data['index'] = new_idx
                    atom_data['chain_idx'] = chain_idx
                    atoms[new_idx] = atom_data
                # Advance offset past the highest atom index in this file
                max_idx = max(parsed.keys())
                offset += max_idx
                chain_idx += 1

    return atoms


def extract_section(filepath, section_name):
    """Extract a section from a GROMACS file"""
    with open(filepath) as f:
        lines = f.readlines()

    in_section = False
    section_lines = []

    for line in lines:
        line = line.rstrip()

        # Check if this is the section header
        if line.strip().startswith(section_name):
            in_section = True
            continue

        # Check if we're at another section
        if line.strip().startswith('[') and in_section:
            break

        # If in section and not a comment, add the line
        if in_section:
            # Remove comments
            comment_pos = line.find(';')
            if comment_pos >= 0:
                line = line[:comment_pos]
            line = line.strip()
            if line:
                section_lines.append(line)

    return section_lines


def parse_atoms_section(section_lines):
    """Parse the [ atoms ] section of a topology"""
    atoms = {}

    for line in section_lines:
        parts = line.split()
        if not parts:
            continue

        try:
            # Basic format: nr type resnr residue atom cgnr charge mass [typeB chargeB massB]
            atom_data = {
                'index': int(parts[0]),
                'type': parts[1],
                'resnr': int(parts[2]),
                'resname': parts[3],
                'atomname': parts[4],
                'charge': float(parts[6]),
                'mass': float(parts[7]),
                'typeB': None,
                'chargeB': None,
                'massB': None
            }

            # Check for dual topology columns (typeB, chargeB, massB)
            if len(parts) >= 11:
                atom_data['typeB'] = parts[8]
                atom_data['chargeB'] = float(parts[9])
                atom_data['massB'] = float(parts[10])

            atoms[atom_data['index']] = atom_data

        except (ValueError, IndexError):
            # Skip malformed lines
            continue

    return atoms


def extract_chain_info(filepath):
    """
    Extract chain information from topology [ molecules ] section.

    Returns dict: {mol_name -> count}
    """
    chain_info = {}

    with open(filepath) as f:
        lines = f.readlines()

    in_molecules = False
    for line in lines:
        line = line.strip()

        if line.startswith('[ molecules ]'):
            in_molecules = True
            continue

        if in_molecules and line.startswith('['):
            break  # End of molecules section

        if in_molecules and line and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mol_name = parts[0]
                    mol_count = int(parts[1])
                    chain_info[mol_name] = mol_count
                except (ValueError, IndexError):
                    pass

    return chain_info


def extract_moleculetype_atoms(itp_file):
    """
    Extract residue ranges for a specific moleculetype from .itp file.

    Returns list of residue numbers in this molecule.
    """
    residues = set()

    with open(itp_file) as f:
        lines = f.readlines()

    in_atoms = False
    for line in lines:
        line_stripped = line.strip()

        if line_stripped.startswith('[ atoms ]'):
            in_atoms = True
            continue

        if in_atoms and line_stripped.startswith('['):
            break

        if in_atoms and line_stripped and not line_stripped.startswith(';'):
            parts = line_stripped.split()
            if len(parts) >= 3:
                try:
                    resnr = int(parts[2])
                    residues.add(resnr)
                except (ValueError, IndexError):
                    pass

    return sorted(residues)
