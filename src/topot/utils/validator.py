"""Strict input validation for the topot CLI.

Each check fails fast with a clear, three-line error message: what failed,
the offending value, and a single concrete Fix line. The CLI catches the
ValidationError, prints it, and exits with status 1.
"""

import re
from pathlib import Path


# Known force-field include basenames that should never be required locally.
# These ship with GROMACS or pmx (forcefield.itp, water/ion models, position
# restraints). The validator skips them when checking local sibling .itp files.
_FF_INCLUDE_BASENAMES = {
    'forcefield.itp', 'tip3p.itp', 'tip4p.itp', 'tip4pew.itp',
    'spc.itp', 'spce.itp', 'ions.itp',
}

# Best-effort atoms-per-residue table for common solvent and ions.
# Used by the atom-count consistency check when the topology has no local
# .itp for a residue (which is the normal case for SOL, NA, CL, etc).
_SOLVENT_ION_ATOMS = {
    'SOL': 3, 'WAT': 3, 'HOH': 3, 'TIP3P': 3, 'TIP3': 3, 'TIP4': 4, 'SPC': 3,
    'NA': 1, 'CL': 1, 'K': 1, 'MG': 1, 'ZN': 1, 'CA': 1, 'FE': 1,
    'ACL': 1, 'ANA': 1,
}


class ValidationError(Exception):
    """Raised by any check when input fails strict validation."""


def _err(what, *detail_lines, fix):
    """Format a three-section error message and raise ValidationError."""
    lines = [what]
    for d in detail_lines:
        lines.append(f"  {d}")
    lines.append(f"  Fix: {fix}")
    raise ValidationError("\n".join(lines))


def _validate_gro_header(gro_file):
    """Check GRO header is well-formed and atom record count matches the
    declared count on line 2. Returns the declared atom count."""
    with open(gro_file) as f:
        lines = f.readlines()

    if len(lines) < 3:
        _err(
            f"ERROR: GRO file too short: {gro_file}",
            f"Found {len(lines)} lines (need title + count + at least one atom)",
            fix="provide a valid GROMACS .gro file with title, atom count, and atom records.",
        )

    try:
        declared = int(lines[1].strip())
    except ValueError:
        _err(
            f"ERROR: GRO atom count line is not an integer: {gro_file}",
            f"Line 2 reads: {lines[1].rstrip()!r}",
            fix="line 2 of a .gro file must be a positive integer (the atom count).",
        )

    if declared <= 0:
        _err(
            f"ERROR: GRO declared atom count is not positive: {gro_file}",
            f"Line 2 declares: {declared}",
            fix="line 2 of a .gro file must be a positive integer (the atom count).",
        )

    # Last line is box vectors; atom records are lines[2 : 2 + declared].
    # The actual record count is bounded by how many lines exist after line 2,
    # excluding the trailing box-vectors line.
    body = lines[2:]
    # Strip trailing blank lines and the box-vectors line (last non-blank).
    while body and not body[-1].strip():
        body.pop()
    if body:
        # The final non-blank line is the box vectors; remove it from records.
        body = body[:-1]

    actual = len(body)
    if actual != declared:
        _err(
            f"ERROR: GRO atom count mismatch: {gro_file}",
            f"Declared on line 2: {declared} atoms",
            f"Found in file:      {actual} atom records",
            fix="ensure the integer on line 2 matches the number of atom records.",
        )
    return declared


def _read_top(top_file):
    """Read the .top file once and return its lines for cheap reuse."""
    with open(top_file) as f:
        return f.read()


_INCLUDE_RE = re.compile(r'#include\s*["\']([^"\']+)["\']')
_MOLTYPE_RE = re.compile(r'^\s*\[\s*moleculetype\s*\]', re.MULTILINE)
_MOLECULES_RE = re.compile(r'^\s*\[\s*molecules\s*\]', re.MULTILINE)


def _extract_molecules_section(content):
    """Return ordered list of (mol_name, count) from [ molecules ]."""
    molecules = []
    in_section = False
    for raw in content.splitlines():
        line = raw.strip()
        if _MOLECULES_RE.match(raw):
            in_section = True
            continue
        if in_section and line.startswith('['):
            break
        if in_section and line and not line.startswith(';'):
            # Strip inline comment.
            line = line.split(';', 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    molecules.append((parts[0], int(parts[1])))
                except ValueError:
                    continue
    return molecules


def _validate_top_structure(top_file):
    """Check the .top contains [ molecules ] and has at least one source
    of [ atoms ] (either #include or local [ moleculetype ])."""
    content = _read_top(top_file)

    if not _MOLECULES_RE.search(content):
        _err(
            f"ERROR: Topology missing [ molecules ] section: {top_file}",
            "TOPOT requires a [ molecules ] section listing every molecule "
            "in the system.",
            fix="append a [ molecules ] section to the .top with each chain, "
                "solvent, and ion entry in the order they appear in the .gro.",
        )

    molecules = _extract_molecules_section(content)
    if not molecules:
        _err(
            f"ERROR: Topology [ molecules ] section is empty: {top_file}",
            "Section exists but contains no parseable 'name count' entries.",
            fix="add lines like 'Protein_chain_A  1' under [ molecules ].",
        )

    includes = _INCLUDE_RE.findall(content)
    has_local_moltype = bool(_MOLTYPE_RE.search(content))

    if not includes and not has_local_moltype:
        _err(
            f"ERROR: Topology has no #include directives and no local [ moleculetype ]: {top_file}",
            "Cannot find any source of [ atoms ] data.",
            fix="either inline a [ moleculetype ] block with [ atoms ], or "
                "#include the per-chain .itp files (e.g. pmx_topol_Protein_chain_A.itp).",
        )

    return {'molecules': molecules, 'includes': includes,
            'has_local_moltype': has_local_moltype}


def _is_ff_include(inc_path):
    """Return True if an #include refers to a force-field file we should
    not require to exist locally."""
    parts = Path(inc_path).parts
    # Any path component ending in .ff is a FF subdirectory reference.
    if any(p.endswith('.ff') for p in parts):
        return True
    # Bare known-FF filenames (e.g. tip3p.itp imported directly).
    if Path(inc_path).name in _FF_INCLUDE_BASENAMES:
        return True
    # Position-restraint files commonly named posre*.itp; skip — these are
    # auto-generated by gmx and not required for dual-topology analysis.
    if Path(inc_path).name.startswith('posre'):
        return True
    return False


def _validate_local_itp_includes(top_file, includes):
    """For each #include that is NOT a force-field reference, assert the
    file exists in the same directory as the .top. Returns the resolved
    list of local .itp paths."""
    top_dir = Path(top_file).parent
    local_paths = []
    missing = []
    for inc in includes:
        if _is_ff_include(inc):
            continue
        inc_path = (top_dir / inc).resolve()
        # Also reject includes that escape top_dir via ../ (treat as FF).
        try:
            inc_path.relative_to(top_dir.resolve())
        except ValueError:
            continue
        if inc_path.exists():
            local_paths.append(inc_path)
        else:
            missing.append((inc, inc_path))

    if missing:
        inc, inc_path = missing[0]
        _err(
            f"ERROR: Missing local #include in {top_file}",
            f"Topology references: {inc}",
            f"Expected location:   {inc_path} (not found)",
            fix="copy the .itp into the same directory as the .top, or fix "
                "the #include path. Force-field files (forcefield.itp, tip3p.itp, "
                "ions.itp, etc.) do not need to be present locally.",
        )
    return local_paths


def _count_atoms_in_itp(itp_path):
    """Return the number of [ atoms ] records in an .itp moleculetype.
    Returns None if no [ atoms ] section is found."""
    with open(itp_path) as f:
        lines = f.readlines()
    in_section = False
    count = 0
    for raw in lines:
        stripped = raw.strip()
        if stripped.startswith('[ atoms ]') or stripped.startswith('[atoms]'):
            in_section = True
            continue
        if in_section and stripped.startswith('['):
            break
        if in_section and stripped and not stripped.startswith(';'):
            count += 1
    return count if in_section else None


def _extract_moltype_name(itp_path):
    """Return the moleculetype name declared in an .itp, or None."""
    with open(itp_path) as f:
        lines = f.readlines()
    in_section = False
    for raw in lines:
        stripped = raw.strip()
        if stripped.startswith('[ moleculetype ]') or stripped.startswith('[moleculetype]'):
            in_section = True
            continue
        if in_section and stripped.startswith('['):
            return None
        if in_section and stripped and not stripped.startswith(';'):
            parts = stripped.split()
            if parts:
                return parts[0]
    return None


def _validate_atom_count_consistency(molecules, local_itps, gro_atom_count):
    """Best-effort: sum atoms across [ molecules ] entries and compare to
    gro_atom_count. ERROR only if every molecule resolved (so the total is
    authoritative) and the totals disagree. Otherwise emit a WARNING."""
    name_to_itp_count = {}
    for itp in local_itps:
        name = _extract_moltype_name(itp)
        atom_count = _count_atoms_in_itp(itp)
        if name and atom_count:
            name_to_itp_count[name] = atom_count

    total = 0
    unresolved = []
    for mol_name, mol_count in molecules:
        if mol_name in name_to_itp_count:
            total += name_to_itp_count[mol_name] * mol_count
        elif mol_name.upper() in _SOLVENT_ION_ATOMS:
            total += _SOLVENT_ION_ATOMS[mol_name.upper()] * mol_count
        else:
            unresolved.append(mol_name)

    if unresolved:
        print(f"  WARNING: Cannot verify total atom count "
              f"(unresolved residues: {', '.join(unresolved)}). "
              f"Skipping atom-count consistency check.")
        return

    if total != gro_atom_count:
        _err(
            "ERROR: Topology and GRO atom counts disagree",
            f"Sum from [ molecules ]: {total} atoms",
            f"Found in GRO file:      {gro_atom_count} atom records",
            fix="ensure the [ molecules ] section lists every chain, solvent, "
                "and ion in the order they appear in the .gro, and that each "
                "chain's .itp [ atoms ] count matches.",
        )


def _validate_ndx(ndx_file):
    """Check the index file has at least one [ name ] header."""
    with open(ndx_file) as f:
        content = f.read()
    if not re.search(r'^\s*\[\s*\S.+?\s*\]', content, re.MULTILINE):
        _err(
            f"ERROR: Index file has no group headers: {ndx_file}",
            "A valid GROMACS .ndx file contains at least one '[ name ]' header.",
            fix="provide a real GROMACS index file (gmx make_ndx output), or omit -n to let TOPOT create a fresh one.",
        )


def _validate_ff_dir(ff_dir):
    """Check the force-field directory exists and contains usable data."""
    ff_dir = Path(ff_dir)
    if not ff_dir.exists():
        _err(
            f"ERROR: Force-field directory not found: {ff_dir}",
            fix="pass an existing directory with --ff-dir, or omit it to use the bundled mutff.",
        )
    if not ff_dir.is_dir():
        _err(
            f"ERROR: --ff-dir is not a directory: {ff_dir}",
            fix="point --ff-dir at a directory (containing *.ff/ subdirectories), not a file.",
        )

    # Acceptable: directory containing one or more *.ff/ children,
    # OR the directory itself is a *.ff/ containing forcefield.itp.
    has_ff_children = any(child.is_dir() and child.name.endswith('.ff')
                          for child in ff_dir.iterdir())
    is_ff_dir = ff_dir.name.endswith('.ff') and (ff_dir / 'forcefield.itp').exists()

    if not (has_ff_children or is_ff_dir):
        _err(
            f"ERROR: --ff-dir contains no force fields: {ff_dir}",
            "Expected one or more *.ff/ subdirectories, or a *.ff/ directory "
            "containing forcefield.itp.",
            fix="use the bundled mutff (omit --ff-dir) or point --ff-dir at a "
                "directory with at least one .ff subdirectory.",
        )


def _validate_dual_topology(topology_atoms):
    """Check that at least one parsed atom carries a state-B column."""
    has_dual = any(atom.get('typeB') for atom in topology_atoms.values())
    if not has_dual:
        _err(
            "ERROR: Not a dual-topology file",
            f"Parsed {len(topology_atoms)} atoms but none have a typeB column.",
            fix="provide a topology produced by pmx (with [ atoms ] rows that "
                "include type/typeB/charge/chargeB/mass/massB columns).",
        )


def validate_inputs_pre_parse(gro_file, top_file, ndx_file):
    """Run all pre-parse checks (1, 2, 3, 4, 5 warning, 7). Returns a dict
    with parsed metadata that the CLI can ignore but is useful for tests."""
    gro_atom_count = _validate_gro_header(gro_file)
    top_info = _validate_top_structure(top_file)
    local_itps = _validate_local_itp_includes(top_file, top_info['includes'])
    _validate_atom_count_consistency(top_info['molecules'], local_itps, gro_atom_count)
    if ndx_file:
        _validate_ndx(ndx_file)
    return {
        'gro_atom_count': gro_atom_count,
        'molecules': top_info['molecules'],
        'local_itps': local_itps,
    }
