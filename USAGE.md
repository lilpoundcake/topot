# TOPOT Usage Guide

## Overview

**topot** is a CLI tool that extracts lambda-specific structures from GROMACS dual topology files created by PMX free energy calculations.

Input: A protein structure with dual topology mutation (from PMX)
Output: Four separate structure files (GRO and PDB) for lambda=0 and lambda=1 states

## Installation

### Prerequisites
- Python 3.8+
- Poetry (for package management)
- micromamba or pip

### With Micromamba (Recommended)

```bash
# Create isolated environment
micromamba create -n topot python=3.10 poetry -c conda-forge
micromamba activate topot

# Install topot
cd /path/to/topot
poetry install

# Verify installation
topot --version
```

### With Standard pip

```bash
python3 -m venv topot_env
source topot_env/bin/activate
pip install poetry
cd /path/to/topot
poetry install
```

## Command Line Usage

```bash
topot -g <GRO_FILE> -p <TOPOLOGY_FILE> -n <INDEX_FILE> -o <OUTPUT_DIR> [--ff-dir <FF_DIR>]
```

### Arguments

| Argument | Short | Long | Description | Required |
|----------|-------|------|-------------|----------|
| GRO file | `-g` | `--gro` | Input coordinate file (.gro format) | Yes |
| Topology | `-p` | `--top` | Input topology file (.top format) | Yes |
| Index | `-n` | `--ndx` | Input index file (.ndx format); if provided, appends new groups; if omitted, creates new file | No |
| Output dir | `-o` | `--output` | Output directory for results | Yes |
| FF directory | | `--ff-dir` | Force field directory (default: ./mutff) | No |
| Version | | `--version` | Show version and exit | No |
| Help | `-h` | `--help` | Show help message | No |

### Examples

**Basic usage (append to existing index):**
```bash
topot -g md_mut.gro -p newtop.top -n index.ndx -o ./results/
```

**Create new index file (no -n flag):**
```bash
topot -g md_mut.gro -p newtop.top -o ./results/
# Creates index.ndx with only lambda_0 and lambda_1 groups
```

**With custom force field directory:**
```bash
topot -g md_mut.gro -p newtop.top -n index.ndx -o ./results/ --ff-dir /path/to/mutff
```

**Test case (H_TRP33TYR mutation):**
```bash
cd tests/H_TRP33TYR
topot -g md_mut.gro -p newtop.top -n index.ndx -o results --ff-dir ../../mutff
```

**Test case without existing index:**
```bash
cd tests/H_TRP33TYR
topot -g md_mut.gro -p newtop.top -o results --ff-dir ../../mutff
# Creates new index.ndx with lambda_0 and lambda_1 groups
```

## Output Directory Safety

The tool **never automatically deletes files** from your output directory. When files already exist, you'll be prompted to choose:

```
⚠ Output directory contains existing item(s)

Choose an option:
  (o) Overwrite existing files only
  (s) Create new subfolder (rename output directory)
  (c) Cancel and exit
```

### Options:
- **(o)** Overwrite: Process normally, overwriting only topot output files while keeping your other files
- **(s)** Subfolder: Create a new subfolder (results_1, results_2, etc.) to keep separate runs organized
- **(c)** Cancel: Exit without making any changes, giving you time to back up or organize files

### Example:
```bash
# First run - creates results/ with output
topot -g md.gro -p newtop.top -o results/

# Second run to same directory
topot -g md.gro -p newtop.top -o results/
# Prompts: choose (o)verwrite, (s)ubfolder, or (c)ancel
# If you choose (s), creates results_1/ instead
```

## Output Files

The tool generates up to 7 files in the output directory:

### GRO Files (All Atoms)

**md_lambda_0_WI.gro** (λ=0, with water/ions)
- All atoms in state A
- Excludes dummy atoms marked with "DUM_" type in state A
- Includes water molecules and ions
- Use for: Full system MD simulations

**md_lambda_0.gro** (λ=0, protein only)
- Protein atoms in state A only
- Excludes water and ions
- Use for: Structural analysis, visualization

**md_lambda_1_WI.gro** (λ=1, with water/ions)
- All atoms in state B
- Excludes dummy atoms marked with "DUM_" type in state B
- Includes water molecules and ions
- Use for: Full system MD simulations

**md_lambda_1.gro** (λ=1, protein only)
- Protein atoms in state B only
- Excludes water and ions
- Use for: Structural analysis, visualization

### PDB Files (Protein Only)

**md_lambda_0.pdb** (λ=0, protein)
- Protein structure in state A
- Standard PDB format with ATOM records
- **Coordinates in Angstroms (Å)** (converted from nm)
- Suitable for protein visualization tools (PyMOL, VMD, etc.)

**md_lambda_1.pdb** (λ=1, protein)
- Protein structure in state B
- Standard PDB format with ATOM records
- **Coordinates in Angstroms (Å)** (converted from nm)
- Suitable for protein visualization tools

### Index File

**index.ndx** (Updated)
- Preserves all existing groups from input
- Adds two new groups:
  - `[ lambda_0 ]` - All atoms in λ=0 state
  - `[ lambda_1 ]` - All atoms in λ=1 state
- Use with: GROMACS analysis tools (gmx energy, gmx trajectory, etc.)

## How It Works

The tool implements a 6-step pipeline:

### Step 1: Validate Input Files
- Checks that GRO, topology, and index files exist
- Reports file locations

### Step 2: Identify Force Field
- Parses topology `#include` directives
- Searches for force field in `--ff-dir`
- Locates `mutres.mtp` mutation definition files
- Supports: AMBER, CHARMM, GROMOS variants

### Step 3: Parse Topology
- Reads `[ atoms ]` sections from all .itp files
- Extracts atom types for both state A and state B
- Identifies dual topology atoms with `typeB` column
- Tracks dummy atoms (marked with "DUM_" prefix)

### Step 4: Parse GRO Coordinates
- Reads atomic positions and velocities
- Matches atoms to topology by (residue_number, atom_name)
- Includes solvent/ion atoms not in topology

### Step 5: Identify Mutations
- Compares atoms between state A and state B
- Counts disappearing atoms (state A only)
- Counts appearing atoms (state B only)
- Reports mutation location and severity

### Step 6: Filter and Output
- **Lambda 0 filtering:** Excludes atoms where `type` starts with "DUM_"
- **Lambda 1 filtering:** Excludes atoms where `typeB` starts with "DUM_"
- Generates GRO files (with and without solvent)
- Generates PDB files (protein only)
- Updates index file with lambda groups

## Understanding Dual Topology

### State A (λ=0) - Original State
Represented by the `type` column in topology. Contains:
- Original atoms with their original atom types
- Atoms that become dummy in state B are marked "DUM_" in `typeB`

### State B (λ=1) - Mutant State
Represented by the `typeB` column in topology. Contains:
- Mutated atoms with new atom types
- Atoms that become dummy in state B have "DUM_" prefix in `typeB`

### Example: W2Y Mutation (Tryptophan → Tyrosine)

Atom CD1 (aromatic ring carbon):
```
State A (λ=0): CD1 with type=CW
State B (λ=1): CD1.gone with type=DUM_CW (dummy)
```
Action: Keep in λ=0 output, remove from λ=1 output

Atom DCD1 (dummy tyrosine aromatic carbon):
```
State A (λ=0): DCD1 with type=DUM_CA (dummy)
State B (λ=1): CD1 with type=CA
```
Action: Remove from λ=0 output, keep in λ=1 output

## Residue Type Classification

The tool classifies residues as:
- **Protein:** Standard amino acids + mutations (ALA, GLY, ..., W2Y, etc.)
- **Solvent:** Water molecules (SOL, WAT, HOH, etc.)
- **Ion:** Inorganic ions (Na+, Cl-, K+, etc.)
- **Unknown:** Unrecognized residues (excluded from PDB)

Classification uses:
1. **Primary:** `residuetypes.dat` from detected force field
2. **Secondary:** Hardcoded standard amino acid names
3. **Tertiary:** Common solvent/ion names
4. **Fallback:** Mark as 'unknown'

## File Format Details

### GRO Format

Fixed-width columnar format (GROMACS coordinate file):
```
Line 1:  Title (arbitrary)
Line 2:  Number of atoms (integer)
Lines 3+: ATOM_RECORD (one per line)
Last:    Box vectors
```

ATOM_RECORD format:
```
Field          Width  Format  Content
residue number 5      int5    1-based residue number
residue name   5      string  e.g., "PHE", "W2Y"
atom name      5      string  e.g., "CA", "CD1"
atom number    5      int5    1-based atom number
x              8      real8.3 Position in nm
y              8      real8.3 Position in nm
z              8      real8.3 Position in nm
vx             8      real8.4 Velocity in nm/ps (optional)
vy             8      real8.4 Velocity in nm/ps (optional)
vz             8      real8.4 Velocity in nm/ps (optional)
```

### Topology Format

GROMACS topology with dual topology (state A and B):
```
[ atoms ]
nr  type  resnr  residue  atom  cgnr  charge  mass  [typeB  chargeB  massB]
1   N3    307    THR      N     1     0.1812  14.01
2   H     307    THR      H     2     0.1934  1.008
...
443 DUM_CA 337   W2Y      DCD1  443   0.0000  3.963  CA  -0.1906  12.01
```

### PDB Format

Standard ATOM records:
```
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C
```

## Troubleshooting

### "Force field not detected"
**Solution:** Use explicit `--ff-dir` parameter pointing to mutff directory
```bash
topot -g md.gro -p newtop.top -n index.ndx -o results --ff-dir /full/path/to/mutff
```

### "0 atoms in PDB output"
**Reason:** Residue type classification failed
**Solution:** Check that residuetypes.dat is in force field directory
```bash
ls /path/to/mutff/residuetypes.dat
```

### Output files are empty
**Reason:** Filtering removed all atoms
**Solution:** Verify topology has proper atom type information
```bash
grep "DUM_" /path/to/newtop.top  # Should find some DUM_ atoms
```

### Index file is missing lambda groups
**Reason:** Tool didn't recognize .ndx format
**Solution:** Verify input index file is valid GROMACS index format

## Performance

- **Processing time:** <2 seconds for 84K atoms
- **Memory:** ~150 MB peak
- **Bottleneck:** GRO file parsing (83K lines)
- **Scaling:** Linear with number of atoms O(n)

## References

- GROMACS File Formats: https://manual.gromacs.org/2026.0/reference-manual/topologies/
- PMX Dual Topology: https://pmx.readthedocs.io/
- Free Energy Calculations: https://pmx.readthedocs.io/
