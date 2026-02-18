#!/usr/bin/env python3
"""
TOPOT - CLI for extracting lambda-specific structures from GROMACS dual topology files
Created by PMX free energy calculations.

Follows a strict 6-step pipeline:
1. Identify force field (AMBER/CHARMM/etc)
2. Identify mutated residues
3. Parse atom states from mutres.mtp
4. Copy non-mutated residues
5. Add/filter mutated residue atoms based on lambda state
6. Generate output files (GRO, PDB, NDX)
"""

import argparse
import sys
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Extract lambda-specific structures from GROMACS dual topology files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  topot -g md.gro -p newtop.top -n index.ndx -o ./results/
  topot -g md.gro -p newtop.top -o ./results/
  topot -g md.gro -p newtop.top -n index.ndx -o ./results/ --ff-dir ./mutff
        """
    )

    parser.add_argument('-g', '--gro', required=True, dest='gro_file',
                        help='Input GRO coordinate file')
    parser.add_argument('-p', '--top', required=True, dest='top_file',
                        help='Input topology file')
    parser.add_argument('-n', '--ndx', required=False, dest='ndx_file', default=None,
                        help='Input index file (optional; if not provided, creates new file)')
    parser.add_argument('-o', '--output', required=True, dest='output_dir',
                        help='Output directory')
    parser.add_argument('--ff-dir', dest='ff_dir', default='./mutff',
                        help='Force field directory (default: ./mutff)')
    parser.add_argument('--version', action='version', version='topot 0.1.0')

    args = parser.parse_args()

    # Convert to absolute paths
    gro_file = Path(args.gro_file).resolve()
    top_file = Path(args.top_file).resolve()
    ndx_file = Path(args.ndx_file).resolve() if args.ndx_file else None
    ff_dir = Path(args.ff_dir).resolve()
    output_dir = Path(args.output_dir)

    # Step 1: Validate input files
    print("Step 1: Validating input files...")
    if not gro_file.exists():
        print(f"ERROR: GRO file not found: {gro_file}")
        sys.exit(1)
    if not top_file.exists():
        print(f"ERROR: Topology file not found: {top_file}")
        sys.exit(1)
    if ndx_file and not ndx_file.exists():
        print(f"ERROR: Index file not found: {ndx_file}")
        sys.exit(1)
    print(f"  ✓ GRO file: {gro_file}")
    print(f"  ✓ Topology file: {top_file}")
    if ndx_file:
        print(f"  ✓ Index file: {ndx_file}")
    else:
        print(f"  ℹ Index file: Will create new file (index.ndx)")

    # Step 2: Identify force field
    print("\nStep 2: Identifying force field...")
    from .utils.ff_detector import detect_force_field
    ff_info = detect_force_field(top_file, ff_dir)
    print(f"  Force field: {ff_info['name']}")
    if ff_info.get('mutres_path'):
        print(f"  Mutation file: {ff_info['mutres_path']}")

    # Step 3: Parse topology and identify atoms
    print("\nStep 3: Parsing topology...")
    from .utils.topology_parser import parse_topology, extract_chain_info, extract_moleculetype_atoms
    topology_atoms = parse_topology(top_file)
    print(f"  Total atoms in topology: {len(topology_atoms)}")
    print(f"  Dual topology atoms: {sum(1 for a in topology_atoms.values() if a.get('typeB'))}")

    # Extract chain information
    chain_molecules = extract_chain_info(top_file)
    if chain_molecules:
        print(f"  Chains found: {list(chain_molecules.keys())}")

    # Step 4: Parse GRO coordinates
    print("\nStep 4: Parsing GRO coordinates...")
    from .utils.gro_parser import parse_gro
    gro_atoms = parse_gro(gro_file)
    print(f"  Total atoms in GRO: {len(gro_atoms)}")

    # Step 5: Identify residue types and mutated residues
    print("\nStep 5: Identifying residue types and mutations...")
    from .utils.residue_classifier import classify_residues, identify_mutations
    residue_types = classify_residues(gro_atoms, ff_dir)
    mutations = identify_mutations(topology_atoms, gro_atoms)
    print(f"  Mutations found: {len(mutations)}")
    for mut in mutations:
        print(f"    - Residue {mut['resnr']}: {mut['atoms_disappear']} disappear, {mut['atoms_appear']} appear")

    # Step 6: Process and generate outputs
    print("\nStep 6: Processing and generating outputs...")
    from .utils.processor import process_dual_topology
    from .utils.mutres_parser import parse_mutres_file

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check for existing files
    existing_files = list(output_dir.glob('*'))
    if existing_files:
        print(f"\n  ⚠ Output directory contains {len(existing_files)} existing item(s)")
        print(f"  Directory: {output_dir}")
        for item in existing_files[:5]:  # Show first 5 items
            print(f"    - {item.name}")
        if len(existing_files) > 5:
            print(f"    ... and {len(existing_files) - 5} more items")

        while True:
            print("\n  Choose an option:")
            print("    (o) Overwrite existing files only")
            print("    (s) Create new subfolder (rename output directory)")
            print("    (c) Cancel and exit")
            choice = input("  Enter choice (o/s/c): ").strip().lower()

            if choice == 'c':
                print("\n  Cancelled by user")
                sys.exit(0)
            elif choice == 's':
                # Create a new subfolder with timestamp/number
                import time
                counter = 1
                while True:
                    new_dir = output_dir.parent / f"{output_dir.name}_{counter}"
                    if not new_dir.exists():
                        output_dir = new_dir
                        output_dir.mkdir(parents=True, exist_ok=True)
                        print(f"  ✓ Created new directory: {output_dir}")
                        break
                    counter += 1
                break
            elif choice == 'o':
                print(f"  ✓ Will overwrite existing files")
                break
            else:
                print("  Invalid choice. Please enter o, s, or c")

    # Parse mutation definitions for atom name mappings
    mutres_path = ff_info.get('mutres_path')
    mutation_maps = {}
    if mutres_path:
        all_mutations = parse_mutres_file(mutres_path)
        mutation_maps = all_mutations
        print(f"  Loaded {len(mutation_maps)} mutation definitions from {Path(mutres_path).name}")

    # Process dual topology and generate outputs
    process_dual_topology(
        topology_atoms=topology_atoms,
        gro_atoms=gro_atoms,
        residue_types=residue_types,
        mutations=mutations,
        ff_dir=ff_dir,
        ndx_file=ndx_file,
        output_dir=output_dir,
        chain_molecules=chain_molecules,
        top_file=top_file,
        mutation_maps=mutation_maps
    )

    print("\n✓ Done!")


if __name__ == '__main__':
    main()
