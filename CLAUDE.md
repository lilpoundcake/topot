# CLAUDE.md - TOPOT Project Guide

This file provides comprehensive guidance to Claude Code when working on the TOPOT repository.

## Project Overview

**TOPOT** (GROMACS Dual Topology Extractor) is a Python CLI tool that processes GROMACS dual topology files created by PMX free energy calculations. It extracts lambda-specific protein structures (λ=0 and λ=1 states) from hybrid topology files.

**Key Use Case:** Free energy calculations for protein mutations require dual topology representations. TOPOT simplifies extracting the separate λ=0 and λ=1 structures for analysis and visualization.

## Repository Structure

```
topot/
├── src/topot/                  # Main package
│   ├── __init__.py            # Package version and metadata
│   ├── cli.py                 # CLI entry point (6-step pipeline)
│   └── utils/                 # Core functionality modules
│       ├── ff_detector.py     # Force field detection from topology #include directives
│       ├── topology_parser.py # GROMACS topology parsing (atoms, types, dual states)
│       ├── gro_parser.py      # GROMACS coordinate file parsing
│       ├── residue_classifier.py  # Classify residues (protein/solvent/ion)
│       ├── mutres_parser.py   # Parse mutation definitions from mutres.mtp
│       └── processor.py       # Main processing: filtering and output generation
├── tests/                     # Test data and cases
│   ├── H_TRP33TYR/           # W2Y mutation (84K atoms, 3 chains)
│   ├── A_ARG155ASH-A_ASP177ASH-A_LYS180ASP/  # Triple mutation test
│   └── test_cli.py           # Basic CLI tests
├── mutff/                    # Force field definitions with 650+ mutations
├── pyproject.toml            # Poetry project configuration
├── setup.py                  # setuptools configuration
├── requirements.txt          # Runtime dependencies
├── requirements-dev.txt      # Development dependencies
├── README.md                 # User documentation (overview, quick start)
├── USAGE.md                  # Comprehensive usage guide and examples
└── CLAUDE.md                 # This file
```

## Core Architecture

### 6-Step Pipeline (cli.py)

1. **Validate** - Check input files exist and are accessible
2. **Detect FF** - Identify force field from topology #include directives (handles both absolute and relative paths)
3. **Parse Topology** - Extract atom types with dual topology info (type, typeB, charge, chargeB)
4. **Parse GRO** - Read atomic coordinates and velocities
5. **Identify Mutations** - Compare atoms between states, count disappearing/appearing atoms
6. **Filter & Output** - Generate lambda-specific GRO/PDB/NDX files

### Data Structures

**Merged Atom Dictionary:**
```python
{
    'index': int,              # Atom serial number
    'atomname': str,           # Atom name (e.g., "CA", "CG", "DCD1")
    'resname': str,            # Residue name (e.g., "ALA", "W2Y")
    'resnr': int,              # Residue number
    'x': float, 'y': float, 'z': float,  # Coordinates in nm
    'vx': float, 'vy': float, 'vz': float,  # Velocities (optional)
    'type': str,               # State A atom type
    'typeB': str or None,      # State B atom type (dual topology)
    'charge': float,           # State A charge
    'chargeB': float or None,  # State B charge
    'mass': float,             # State A mass
    'massB': float or None,    # State B mass
}
```

### Filtering Rules

**Lambda 0 (State A):**
- ✓ Include atoms where `type` does NOT start with "DUM_"
- ✗ Exclude dummy atoms (placeholders for atoms appearing in state B)

**Lambda 1 (State B):**
- ✓ Include atoms where `typeB` does NOT start with "DUM_"
- ✗ Exclude dummy atoms (placeholders for atoms disappearing from state A)

### Atom Renaming (for lambda_1)

Uses `mutres.mtp` mutation definition file to rename atoms in lambda_1 state.

Example W2Y (Tryptophan → Tyrosine):
- State A (lambda_0): DCD1 (dummy carbon), HV1 (dummy hydrogen) - shown as is
- State B (lambda_1): Renamed to CD1, HD1 (showing tyrosine aromatic positions)

**Special Handling for D-Prefixed Atoms:**
- Convention: "D" prefix indicates dual topology atom only in state B
- Example: DHD2 → HD2 (strips D prefix and checks if result is valid state B name)
- Applies both direct mappings and fallback D-prefix stripping

## Key Implementation Details

### processor.py Critical Functions

**`build_resnr_to_mutation_map(gro_atoms, mutation_maps)`**
- Creates mapping: residue_number → mutation_name (e.g., 337 → "W2Y")
- Used to apply correct mutation-specific atom renamings

**`filter_lambda_0(atoms)` and `filter_lambda_1(atoms)`**
- Lambda_0: excludes atoms with type starting with "DUM_"
- Lambda_1: excludes atoms with typeB starting with "DUM_" (or uses type if typeB is None)

**`write_pdb(atoms, filepath, ..., resnr_to_mutation=None)`**
- Uses BioPython PDBIO for proper PDB format compliance
- Applies atom renamings ONLY from the specific residue's mutation
- Converts coordinates from nm (GRO) to Å (PDB): multiply by 10.0
- Handles element inference from renamed atom names

**`write_gro(atoms, filepath, title)`**
- Outputs coordinates in GRO format (nm)
- Preserves fixed-width column alignment

**`update_ndx(lambda_0_atoms, lambda_1_atoms, input_ndx, output_ndx, residue_types)`**
- Appends 4 groups to index file:
  - `lambda_0` - All lambda_0 atoms
  - `lambda_1` - All lambda_1 atoms
  - `lambda_0_wo_water_and_ions` - Protein only (λ=0)
  - `lambda_1_wo_water_and_ions` - Protein only (λ=1)
- Preserves all existing groups when appending

### ff_detector.py - Force Field Detection

**Dual Input Mode Support:**
- Can accept `ff_dir` as parent directory (containing .ff subdirectories)
- Can also accept `ff_dir` as a .ff directory itself
- Handles both: `./mutff` and `./mutff/amber99sb-star-ildn-mut.ff`
- Properly extracts FF name from relative paths with `../` notation

### mutres_parser.py - Mutation Definition Parsing

**Parses mutres.mtp files with 650+ mutation definitions**
- Format: `[ mutation_name ]` header followed by `[ morphes ]` section
- Morphes format: `state_A_name state_A_type -> state_B_name state_B_type`
- Stores bidirectional mappings: `atom_A_to_B` and `atom_B_to_A`
- Filters out ".gone" atoms (virtual documentation atoms)

### residue_classifier.py - Residue Classification

**Classifies residues as:**
- Protein: Standard amino acids + mutations (ALA, GLY, ..., W2Y, etc.)
- Solvent: Water molecules (SOL, WAT, HOH, TIP3P, etc.)
- Ion: Inorganic ions (Na+, Cl-, K+, Mg2+, etc.)
- Unknown: Unrecognized residues

**Classification priority:**
1. Primary: `residuetypes.dat` from detected force field (normalized to lowercase)
2. Secondary: Hardcoded standard amino acid names
3. Tertiary: Common solvent/ion names
4. Fallback: Mark as 'unknown'

## Output Directory Safety

The tool **never automatically deletes files** in the output directory. When files exist, users are prompted:
- **(o)** Overwrite: Write output files to existing directory (non-conflicting files preserved)
- **(s)** Subfolder: Create numbered subfolder (results_1, results_2, etc.)
- **(c)** Cancel: Exit without making any changes

**Implementation:** `cli.py` lines 112-145

## Index File Handling

The `-n` (index) flag is now optional:

- **With `-n index.ndx`**: Reads existing index file and appends `lambda_0`, `lambda_1`, `lambda_0_wo_water_and_ions`, `lambda_1_wo_water_and_ions` groups
- **Without `-n`**: Creates new `index.ndx` file with only these 4 groups

Both modes preserve all existing groups when appending.

## Development Workflow

### Running Tests

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest --cov=src tests/

# Run specific test
pytest tests/test_cli.py -v
```

### Code Quality

```bash
# Format code
black src/ tests/

# Check imports
isort src/ tests/

# Lint
flake8 src/ tests/

# Type check
mypy src/
```

### Build Distribution

```bash
# Create wheel and source distribution
poetry build

# Verify wheel
tar -tzf dist/*.tar.gz | head

# Upload to PyPI (when ready)
twine upload dist/*
```

## Common Development Tasks

### Adding a New Parser Module

```python
# src/topot/utils/new_parser.py
"""Parse new file format"""

def parse_newformat(filepath):
    """Parse file and return data structure"""
    data = []
    with open(filepath) as f:
        for line in f:
            # Parse logic
            pass
    return data
```

Then import in `cli.py` and add to the pipeline.

### Modifying Filtering Logic

Edit `filter_lambda_0()` and `filter_lambda_1()` in `processor.py`. Currently:
- Lambda_0: Exclude atoms with `type` starting with "DUM_"
- Lambda_1: Exclude atoms with `typeB` starting with "DUM_"

To change: modify the `if atom_type.startswith('DUM_')` condition.

### Extending for New Output Formats

Add new function in `processor.py`:

```python
def write_newformat(atoms, filepath, chain_map=None):
    """Write atoms in new format"""
    with open(filepath, 'w') as f:
        for atom in atoms:
            # Format-specific writing
            pass
```

Then call from `process_dual_topology()`.

### Testing a Specific Case

```bash
cd tests/H_TRP33TYR

# Run tool with existing index file
topot -g md_mut.gro -p newtop.top -n index.ndx -o results \
      --ff-dir ../../mutff

# Or create new index file
topot -g md_mut.gro -p newtop.top -o results \
      --ff-dir ../../mutff

# Verify outputs
ls results/md_lambda*.{gro,pdb}

# Check atom counts
echo "Lambda_0 protein atoms:"
grep "^ " results/md_lambda_0.gro | wc -l

echo "Lambda_1 protein atoms:"
grep "^ " results/md_lambda_1.gro | wc -l
```

Expected: λ_0 = 6487 atoms, λ_1 = 6484 atoms (3 atom difference for W→Y)

## Dependencies and Versions

**Runtime:**
- Python 3.8+ (tested on 3.8, 3.9, 3.10, 3.11, 3.12)
- numpy >= 1.19
- biopython >= 1.81

**Development:**
- pytest, pytest-cov
- black, flake8, mypy, isort
- build, twine, wheel
- poetry

**No C extensions** - Pure Python for maximum portability.

## Important Notes for Modifications

1. **Atom Name Handling:** Always preserve original atom names from GRO when not applying mutation renames
2. **Coordinate Units:** GRO files use nanometers (nm), PDB uses Ångströms (Å). Convert: multiply by 10.0
3. **Residue Number Mapping:** Some mutations affect multiple atoms; always check `resnr` when applying mutations
4. **Chain Separation:** Uses residue number jumps or topology chain info for detection
5. **Type Vs. TypeB:** In topology, "type" is state A, "typeB" is state B. This is fundamental to dual topology
6. **Element Inference:** PDB element symbols derived from first uppercase letter of atom name (after renaming)

## Known Limitations

1. **Velocity Handling:** Velocities are preserved from GRO but not validated
2. **Box Vectors:** Output GRO box vectors are set to zeros (not propagated from input)
3. **Multiple Mutations:** Assumes max 1 mutation per residue (but handles multi-residue mutations)
4. **Force Field Parsing:** Only detects FF from #include directives, not from GMXLIB environment variable

## Future Enhancement Ideas

1. Add support for trajectory files (XTC format)
2. Preserve box vectors from input GRO
3. Add validation mode (check consistency between topology and GRO)
4. Add performance profiling and optimization
5. Add support for gmx convert-tpr intermediate format
6. Add batch processing for multiple mutations

## Quick Reference

### File Paths to Know

- **CLI Entry:** `src/topot/cli.py`
- **Main Processing:** `src/topot/utils/processor.py`
- **Force Field Detection:** `src/topot/utils/ff_detector.py`
- **Test Case 1:** `tests/H_TRP33TYR/` (W2Y single mutation)
- **Test Case 2:** `tests/A_ARG155ASH-A_ASP177ASH-A_LYS180ASP/` (Triple mutation)
- **Force Fields:** `mutff/` (650+ mutations across AMBER, CHARMM, GROMOS variants)

### Key Classes/Functions

- `parse_topology()` - Reads GROMACS topology files
- `parse_gro()` - Reads GROMACS coordinate files
- `merge_atoms()` - Combines topology and coordinate data
- `classify_residues()` - Classifies residues by type
- `process_dual_topology()` - Main processing function
- `filter_lambda_0/1()` - Filters atoms for each state
- `write_pdb/gro()` - Output file generation
- `update_ndx()` - Updates index files

### Error Messages

- "Force field not detected" → Pass explicit `--ff-dir`
- "0 atoms in PDB" → Check residuetypes.dat exists in FF directory
- "Empty output" → Verify topology has DUM_ atoms for filtering
- "Index group mismatch" → Ensure input .ndx file is valid GROMACS format

## References

- [GROMACS File Formats](https://manual.gromacs.org/2026.0/reference-manual/topologies/)
- [PMX Documentation](https://pmx.readthedocs.io/)
- [Dual Topology FE Calculations](https://pmx.readthedocs.io/guide/hybrid-topology.html)
- [BioPython PDBIO](https://biopython.org/)
- [GROMACS Topology Manual](https://manual.gromacs.org/current/reference-manual/file-formats.html#top)

## Current Status

✅ Core functionality complete
✅ Single and multi-mutation support
✅ PDB output with proper coordinates and element symbols
✅ Protein-only index groups (wo_water_and_ions)
✅ Output directory safety (user-controlled overwrite/subfolder/cancel)
✅ Optional index file handling (create new or append)
✅ Smart atom renaming with D-prefix fallback
✅ All tests passing
✅ Wheel distribution created
✅ Full documentation (README.md, USAGE.md, CLAUDE.md)
