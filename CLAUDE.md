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
│   ├── data/                  # Bundled data files (included in wheel)
│   │   └── mutff/             # Force field definitions with 650+ mutations
│   │       ├── amber14sbmut.ff/
│   │       ├── amber99sb-star-ildn-mut.ff/
│   │       ├── amber99sb-star-ildn-bsc1-mut.ff/
│   │       ├── amber99sb-star-ildn-dna-mut.ff/
│   │       ├── charmm22star-mut.ff/
│   │       ├── charmm36m-mut.ff/
│   │       └── residuetypes.dat, atomtypes.atp, etc.
│   └── utils/                 # Core functionality modules
│       ├── ff_detector.py     # Force field detection from topology #include directives
│       ├── topology_parser.py # GROMACS topology parsing (atoms, types, dual states)
│       ├── gro_parser.py      # GROMACS coordinate file parsing
│       ├── residue_classifier.py  # Classify residues (protein/solvent/ion)
│       ├── mutres_parser.py   # Parse mutation definitions from mutres.mtp
│       └── processor.py       # Main processing: filtering and output generation
├── tests/                     # Test data and cases
│   ├── H_TRP33TYR/           # W2Y mutation (84K atoms, 3 chains)
│   ├── L_ASN57HID-H_TYR104GLN/  # N2H + Y2Q dual mutation across chains (34K atoms)
│   └── A_ARG155ASH-A_ASP177ASH-A_LYS180ASP/  # Triple mutation test
├── pyproject.toml            # Poetry project configuration (includes data files)
├── setup.py                  # setuptools configuration (includes data files)
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

### Topology Parsing - Multi-Chain .itp Handling

Each `.itp` file may restart atom numbering from 1. `parse_topology()` applies an offset
to atom indices when combining atoms from multiple files, ensuring no chain's data is lost.

**Critical:** Without the offset, `atoms.update(parsed)` overwrites previous chains' atoms,
causing all dual topology info (type/typeB) to be lost for those chains. This was fixed to
handle arbitrary numbers of chains with independent atom numbering.

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

### ff_detector.py - Force Field Detection (Multi-Strategy)

**Three-tier detection strategy:**

1. **Primary Method: #include Directive Parsing**
   - Reads `#include` statements from topology file
   - Extracts FF name from paths (e.g., `#include "amber99sb-star-ildn-mut.ff/forcefield.itp"`)
   - Handles both absolute and relative paths with `../` notation
   - Direct match: fastest and most reliable

2. **Fallback 1: Similar Name Matching**
   - If exact FF not found, tries stripping common suffixes
   - Examples:
     - `amber99sb-star-ildn-mut-charged` → tries `amber99sb-star-ildn-mut`
     - `amber99sb-star-ildn-mut-dna` → tries `amber99sb-star-ildn-mut`
   - Handles cases where topology specifies variant that doesn't exist in `./mutff`

3. **Fallback 2: Atom Type Inference**
   - Extracts atom types from `[ atoms ]` section in topology
   - Parses `atomtypes.atp` from each available FF
   - Calculates overlap: must have >50% atom types in common
   - Used when topology has no `#include` directives

**Dual Input Mode Support:**
- Can accept `ff_dir` as parent directory (containing .ff subdirectories)
- Can also accept `ff_dir` as a .ff directory itself
- Handles both: `./mutff` and `./mutff/amber99sb-star-ildn-mut.ff`

**Return Dictionary:**
```python
{
    'name': str,                    # Force field name (or 'unknown')
    'path': str,                    # Full path to FF directory
    'mutres_path': str or None,     # Path to mutres.mtp if exists
    'detection_method': str,        # How FF was detected
    'available_ff': list            # FFs available in ff_dir
}
```

### Bundled Force Field Data

**Distribution Strategy:**
- Force field data (6 variants, 650+ mutations) is bundled inside the Python wheel
- Located at `src/topot/data/mutff/` during development
- After `pip install topot`, accessible via `Path(__file__).parent / 'data' / 'mutff'`
- Makes tool self-contained: no need for `--ff-dir` after installation

**Default Resolution (cli.py):**
```python
_BUNDLED_FF_DIR = Path(__file__).parent / 'data' / 'mutff'
parser.add_argument('--ff-dir', dest='ff_dir', default=str(_BUNDLED_FF_DIR), ...)
```
- `Path(__file__)` always points to installed package location
- Works correctly after `pip install` (even if git repo is deleted)
- `--ff-dir` override still works for custom force fields

**Wheel Configuration:**
- `pyproject.toml`: adds `include` directive for Poetry
- `setup.py`: adds `package_data` and `include_package_data` for setuptools
- Wheel size: 2.0 MB (compressed), includes all 217 FF data files

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

# Run tool with existing index file (mutff is bundled, no --ff-dir needed)
topot -g md_mut.gro -p newtop.top -n index.ndx -o results

# Or create new index file
topot -g md_mut.gro -p newtop.top -o results

# Or with custom FF directory override
topot -g md_mut.gro -p newtop.top -o results --ff-dir /custom/path/mutff

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
7. **Multi-Chain .itp Files:** Each `.itp` restarts atom numbering from 1; `parse_topology()` offsets indices to avoid collisions

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

- **CLI Entry:** `src/topot/cli.py` (default FF dir resolved here)
- **Main Processing:** `src/topot/utils/processor.py`
- **Force Field Detection:** `src/topot/utils/ff_detector.py`
- **Test Case 1:** `tests/H_TRP33TYR/` (W2Y single mutation, 84K atoms)
- **Test Case 2:** `tests/L_ASN57HID-H_TYR104GLN/` (N2H+Y2Q dual mutation across chains, 34K atoms)
- **Test Case 3:** `tests/A_ARG155ASH-A_ASP177ASH-A_LYS180ASP/` (Triple mutation)
- **Bundled Force Fields:** `src/topot/data/mutff/` (650+ mutations, included in wheel)

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

- "Force field not detected" → Check topology #include directives, or pass explicit `--ff-dir`
- "0 atoms in PDB" → Verify residuetypes.dat exists (bundled or custom --ff-dir)
- "Empty output" → Verify topology has DUM_ atoms for filtering
- "Index group mismatch" → Ensure input .ndx file is valid GROMACS format
- "Mutation file not found" → Check mutres.mtp exists in detected FF (bundled or --ff-dir)
- "Atom X defined twice in residue" → Topology parsing lost chain data; check .itp atom index offsets

## References

- [GROMACS File Formats](https://manual.gromacs.org/2026.0/reference-manual/topologies/)
- [PMX Documentation](https://pmx.readthedocs.io/)
- [Dual Topology FE Calculations](https://pmx.readthedocs.io/guide/hybrid-topology.html)
- [BioPython PDBIO](https://biopython.org/)
- [GROMACS Topology Manual](https://manual.gromacs.org/current/reference-manual/file-formats.html#top)

## Current Status

✅ Core functionality complete
✅ Single and multi-mutation support
✅ Multi-chain topology parsing (atom index offset for per-chain .itp files)
✅ PDB output with proper coordinates and element symbols
✅ Protein-only index groups (wo_water_and_ions)
✅ Output directory safety (user-controlled overwrite/subfolder/cancel)
✅ Optional index file handling (create new or append)
✅ Smart atom renaming with D-prefix fallback
✅ Three-tier force field detection (#include → similar names → atom type inference)
✅ Bundled force field data (6 FF variants, 650+ mutations) in wheel
✅ Optional `-o` output directory (defaults to `<gro_stem>_topot`)
✅ All test cases passing (3 test cases: single, dual cross-chain, triple mutation)
✅ Wheel distribution created (2.0 MB with bundled data)
✅ Full documentation (README.md, USAGE.md, CLAUDE.md)
