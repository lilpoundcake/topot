# TOPOT - GROMACS Dual Topology Extractor

A CLI tool for extracting lambda-specific structures from GROMACS dual topology files created by PMX free energy calculations.

## Quick Start

```bash
# Install in isolated environment
micromamba create -n topot python=3.10 poetry -c conda-forge
micromamba activate topot

cd /path/to/topot
poetry install

# Run (mutff is bundled, no --ff-dir needed)
topot -g md_mut.gro -p newtop.top -n index.ndx -o ./results/
```

## What It Does

Takes a protein structure with dual topology mutation (from PMX) and generates:
- 4 GRO files: λ=0 and λ=1, with and without solvent
- 2 PDB files: λ=0 and λ=1 (protein only)
- 1 NDX file: Updated with lambda groups

## Example: W2Y Mutation

Input: Tryptophan (W) mutating to Tyrosine (Y) at position 337
```
λ=0 Output: All tryptophan atoms + backbone (DUM tyrosine atoms removed)
λ=1 Output: All tyrosine atoms + backbone (DUM tryptophan atoms removed)
```

## Key Features

✅ **Bundled force field data** (6 FF variants with 650+ mutations)
✅ **Automatic force field detection** (AMBER, CHARMM, GROMOS)
✅ **Works with ANY mutation type** (general-purpose solution)
✅ **Intelligent atom filtering** based on dual topology states
✅ **Preserves atom connectivity** (no reordering)
✅ **Multiple output formats** (GRO + PDB + NDX)
✅ **Sub-2-second processing** for 84K atom systems
✅ **Zero additional dependencies** (pure Python + numpy/biopython)

## Installation

### With Micromamba (Recommended)
```bash
micromamba create -n topot python=3.10 poetry -c conda-forge
micromamba activate topot
cd /path/to/topot
poetry install
```

### With Pip
```bash
python3 -m venv topot_env
source topot_env/bin/activate
pip install poetry
cd /path/to/topot
poetry install
```

## Usage

```bash
topot -g <GRO_FILE> -p <TOPOLOGY_FILE> -n <INDEX_FILE> -o <OUTPUT_DIR> [--ff-dir <FF_DIR>]
```

### Arguments
- `-g, --gro FILE` - Input GROMACS coordinate file (required)
- `-p, --top FILE` - Input topology file (required)
- `-n, --ndx FILE` - Input index file (optional; if provided, appends to it; if not provided, creates new)
- `-o, --output DIR` - Output directory (optional; default: `<gro_stem>_topot` in current directory)
- `--ff-dir DIR` - Force field directory (optional; default: bundled mutff with AMBER/CHARMM/GROMOS)
- `--version` - Show version
- `-h, --help` - Show help

### Example
```bash
# Test case: W2Y mutation in protein
cd tests/H_TRP33TYR
topot -g md_mut.gro -p newtop.top -n index.ndx -o results
# (mutff is bundled, no --ff-dir needed)

# Or with explicit force field directory
topot -g md_mut.gro -p newtop.top -n index.ndx -o results --ff-dir /custom/path/mutff
```

## Output Files

Generated in output directory:

| File | Contains | Use |
|------|----------|-----|
| `md_lambda_0_WI.gro` | All atoms at λ=0 | Full system with solvent |
| `md_lambda_0.gro` | Protein at λ=0 | Protein structure analysis |
| `md_lambda_1_WI.gro` | All atoms at λ=1 | Full system with solvent |
| `md_lambda_1.gro` | Protein at λ=1 | Protein structure analysis |
| `md_lambda_0.pdb` | Protein at λ=0 | Visualization in PyMOL/VMD |
| `md_lambda_1.pdb` | Protein at λ=1 | Visualization in PyMOL/VMD |
| `index.ndx` | Updated groups | GROMACS analysis tools |

## Output Directory Handling

The tool **never automatically deletes files** in the output directory. When files already exist:

1. **Overwrite (o)** - Process and write output files while preserving non-conflicting files
2. **Create subfolder (s)** - Create a numbered subfolder (results_1, results_2, etc.) instead
3. **Cancel (c)** - Exit without making any changes

This ensures your existing analysis and notes are always preserved.

## How It Works

### 6-Step Pipeline

1. **Validate** - Check input files exist
2. **Detect FF** - Identify force field from topology #include directives
3. **Parse Topology** - Read atom types and dual topology (state A and B)
4. **Parse Coordinates** - Read atomic positions and velocities
5. **Identify Mutations** - Compare atoms between states, count changes
6. **Filter & Output** - Generate lambda-specific files

### Filtering Rules

**Lambda 0 (State A):**
- ✓ Include: atoms with regular type
- ✗ Exclude: atoms with type starting with "DUM_"

**Lambda 1 (State B):**
- ✓ Include: atoms with regular typeB
- ✗ Exclude: atoms with typeB starting with "DUM_"

### Example Filtering

Tryptophan → Tyrosine (W2Y):
```
CD1 (Trp aromatic)
  State A: CD1 with type=CW
  State B: CD1.gone with type=DUM_CW (dummy)
  → Keep in λ=0, Remove in λ=1

DCD1 (Tyr aromatic)
  State A: DCD1 with type=DUM_CA (dummy)
  State B: CD1 with type=CA
  → Remove in λ=0, Keep in λ=1
```

## File Formats

### Input: GRO Format
GROMACS coordinate file with fixed-width columns:
```
Title
N_atoms
 resnr resname atomname atomid x y z [vx vy vz]
...
0 0 0 0
```

### Input: Topology Format
GROMACS topology with dual states:
```
[ atoms ]
nr  type  resnr  residue  atom  cgnr  charge  mass  [typeB  chargeB  massB]
1   N3    307    THR      N     1     0.181   14.0
```

### Output: GRO Format
Same as input, filtered for lambda state

### Output: PDB Format
Standard ATOM records:
```
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C
```

### Output: Index Format
GROMACS index file with groups:
```
[ lambda_0 ]
1 2 3 4 5 ...

[ lambda_1 ]
1 2 3 4 5 ...
```

## Performance

- **Processing time:** <2 seconds for 84,949 atom structure
- **Memory:** ~150 MB peak
- **Bottleneck:** GRO file parsing (fixed-width format)
- **Scaling:** Linear with atom count O(n)

## Troubleshooting

### "Force field not detected"
Use explicit FF directory:
```bash
topot ... --ff-dir /full/path/to/mutff
```

### "0 atoms in output"
Check residue classification:
```bash
# Verify residuetypes.dat exists
ls /path/to/mutff/residuetypes.dat
```

### "Empty protein PDB"
Verify topology has proper atom type information:
```bash
grep "DUM_" /path/to/newtop.top  # Should find dummy atoms
```

## Dependencies

**Runtime:**
- Python 3.8+
- numpy >= 1.19
- biopython >= 1.81

**Development:**
- pytest (testing)
- pytest-cov (coverage)
- black (code formatting)
- flake8 (linting)
- mypy (type checking)
- isort (import sorting)

## Documentation

- **[README.md](README.md)** - Overview and quick start (this file)
- **[USAGE.md](USAGE.md)** - Complete usage guide with examples
- **[CLAUDE.md](CLAUDE.md)** - Technical implementation guide for developers

## Project Structure

```
topot/
├── src/topot/              # Main package
│   ├── __init__.py
│   ├── cli.py              # CLI entry point
│   ├── data/mutff/         # Bundled force field data (6 FF + 650+ mutations)
│   │   ├── amber14sbmut.ff/
│   │   ├── amber99sb-star-ildn-mut.ff/
│   │   ├── amber99sb-star-ildn-bsc1-mut.ff/
│   │   ├── amber99sb-star-ildn-dna-mut.ff/
│   │   ├── charmm22star-mut.ff/
│   │   ├── charmm36m-mut.ff/
│   │   ├── residuetypes.dat
│   │   └── ... (atomtypes, elements, etc.)
│   └── utils/              # Core functionality
│       ├── ff_detector.py
│       ├── topology_parser.py
│       ├── gro_parser.py
│       ├── residue_classifier.py
│       └── processor.py
├── tests/                  # Test data and cases
│   ├── H_TRP33TYR/         # W2Y single mutation (84K atoms, 3 chains)
│   ├── L_ASN57HID-H_TYR104GLN/  # Dual mutation across chains (34K atoms)
│   └── A_ARG155ASH-A_ASP177ASH-A_LYS180ASP/  # Triple mutation
├── pyproject.toml          # Poetry configuration (includes data files)
├── setup.py                # setuptools configuration
├── README.md               # This file
├── USAGE.md                # Comprehensive usage guide
└── CLAUDE.md               # Technical implementation guide
```

## Test Cases

### H_TRP33TYR (Single mutation, 84K atoms)
Chain H, Tryptophan → Tyrosine at position 337

```bash
cd tests/H_TRP33TYR
topot -g md_mut.gro -p topol.top -o results
```

Expected: λ=0: 6,487 protein atoms, λ=1: 6,484 protein atoms

### L_ASN57HID-H_TYR104GLN (Dual mutation across chains, 34K atoms)
Chain L, Asparagine → Histidine at position 57; Chain H, Tyrosine → Glutamine at position 215

```bash
cd tests/L_ASN57HID-H_TYR104GLN
topot -g md.gro -p topol.top -o results
```

Expected: λ=0: 3,802 protein atoms, λ=1: 3,801 protein atoms

## References

- [GROMACS File Formats](https://manual.gromacs.org/2026.0/reference-manual/topologies/)
- [PMX Documentation](https://pmx.readthedocs.io/)
- [Dual Topology FE Calculations](https://pmx.readthedocs.io/)

## License

[Your License Here]

## Contributing

Contributions welcome! Please:
1. Follow PEP 8 style guide
2. Add type hints to all functions
3. Include docstrings
4. Run tests: `pytest tests/`
5. Format code: `black src/`
6. Type check: `mypy src/`

## Support

For issues or questions:
1. Check [USAGE.md](USAGE.md) for common problems
2. Review [CLAUDE.md](CLAUDE.md) for technical details
3. Check test cases in `tests/H_TRP33TYR/`, `tests/L_ASN57HID-H_TYR104GLN/`, or `tests/A_ARG155ASH-A_ASP177ASH-A_LYS180ASP/`
