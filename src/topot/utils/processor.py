"""Main processing logic for dual topology extraction"""

from pathlib import Path
import numpy as np
from Bio.PDB import Structure, Model, Chain, Residue, Atom, PDBIO


def build_resnr_to_mutation_map(gro_atoms, mutation_maps):
    """
    Build a mapping from residue number to mutation name.

    Uses the resname from GRO atoms to determine which mutation each residue has.
    Only maps residues that have a corresponding entry in mutation_maps.

    Returns dict: {resnr -> mutation_name}
    """
    resnr_to_mutation = {}

    if not mutation_maps:
        return resnr_to_mutation

    # Create a set of known mutation names for quick lookup
    known_mutations = set(mutation_maps.keys())

    # Map residue numbers to their resnames
    resnr_to_resname = {}
    for atom in gro_atoms:
        resnr = atom['resnr']
        if resnr not in resnr_to_resname:
            resnr_to_resname[resnr] = atom['resname']

    # Map residue numbers to mutations based on resname
    for resnr, resname in resnr_to_resname.items():
        if resname in known_mutations:
            resnr_to_mutation[resnr] = resname

    return resnr_to_mutation


def process_dual_topology(topology_atoms, gro_atoms, residue_types, mutations,
                          ff_dir, ndx_file, output_dir, chain_molecules=None, top_file=None,
                          mutation_maps=None):
    """
    Main processing pipeline:
    1. Filter atoms for lambda_0 and lambda_1
    2. Write GRO files (with and without solvent)
    3. Write PDB files (protein only)
    4. Update NDX file
    """

    # Build residue number to mutation name mapping
    resnr_to_mutation = build_resnr_to_mutation_map(gro_atoms, mutation_maps)

    # Merge topology and GRO atoms
    merged_atoms = merge_atoms(topology_atoms, gro_atoms)

    # Filter atoms for each lambda state
    lambda_0_atoms = filter_lambda_0(merged_atoms)
    lambda_1_atoms = filter_lambda_1(merged_atoms)

    # Separate by residue type (case-insensitive)
    lambda_0_with_solvent = lambda_0_atoms
    lambda_0_protein_only = [a for a in lambda_0_atoms
                              if residue_types.get(a.get('resnr'), 'unknown').lower() == 'protein']

    lambda_1_with_solvent = lambda_1_atoms
    lambda_1_protein_only = [a for a in lambda_1_atoms
                              if residue_types.get(a.get('resnr'), 'unknown').lower() == 'protein']

    # Write GRO files
    # _WI versions include water and ions, regular versions do not
    print(f"  Writing GRO files...")
    write_gro(lambda_0_with_solvent, output_dir / 'md_lambda_0_WI.gro', 'PMX MODEL - Lambda 0 with water and ions')
    write_gro(lambda_0_protein_only, output_dir / 'md_lambda_0.gro', 'PMX MODEL - Lambda 0')
    write_gro(lambda_1_with_solvent, output_dir / 'md_lambda_1_WI.gro', 'PMX MODEL - Lambda 1 with water and ions')
    write_gro(lambda_1_protein_only, output_dir / 'md_lambda_1.gro', 'PMX MODEL - Lambda 1')

    print(f"    md_lambda_0_WI.gro: {len(lambda_0_with_solvent)} atoms (all)")
    print(f"    md_lambda_0.gro: {len(lambda_0_protein_only)} atoms (protein only)")
    print(f"    md_lambda_1_WI.gro: {len(lambda_1_with_solvent)} atoms (all)")
    print(f"    md_lambda_1.gro: {len(lambda_1_protein_only)} atoms (protein only)")

    # Build chain map using actual topology information
    chain_map = build_chain_map(lambda_0_protein_only, chain_molecules=chain_molecules,
                                topology_atoms=topology_atoms, top_file=top_file)

    # Write PDB files
    print(f"  Writing PDB files...")
    write_pdb(lambda_0_protein_only, output_dir / 'md_lambda_0.pdb', 'Lambda 0 - Protein',
              chain_map=chain_map, mutation_maps=mutation_maps, lambda_state=0,
              topology_atoms=topology_atoms, resnr_to_mutation=resnr_to_mutation)
    write_pdb(lambda_1_protein_only, output_dir / 'md_lambda_1.pdb', 'Lambda 1 - Protein',
              chain_map=chain_map, mutation_maps=mutation_maps, lambda_state=1,
              topology_atoms=topology_atoms, resnr_to_mutation=resnr_to_mutation)

    print(f"    md_lambda_0.pdb: {len(lambda_0_protein_only)} atoms")
    print(f"    md_lambda_1.pdb: {len(lambda_1_protein_only)} atoms")

    # Update or create NDX file
    print(f"  Updating NDX file...")
    update_ndx(lambda_0_atoms, lambda_1_atoms, ndx_file, output_dir / 'index.ndx', residue_types)
    if ndx_file:
        print(f"    Added lambda_0, lambda_1, and protein-only groups to existing index")
    else:
        print(f"    Created new index file with lambda_0, lambda_1, and protein-only groups")


def merge_atoms(topology_atoms, gro_atoms):
    """
    Merge topology information with GRO coordinates.

    Returns list of atoms with combined information.
    """
    merged = []

    # Create a mapping from (resnr, atomname) to topology info
    topology_map = {}
    for idx, topo_atom in topology_atoms.items():
        key = (topo_atom['resnr'], topo_atom['atomname'])
        topology_map[key] = topo_atom

    # Match GRO atoms with topology
    for gro_atom in gro_atoms:
        key = (gro_atom['resnr'], gro_atom['atomname'])

        atom = {**gro_atom}  # Start with GRO data

        # Add topology information if available
        if key in topology_map:
            topo = topology_map[key]
            atom['type'] = topo['type']
            atom['typeB'] = topo['typeB']
            atom['charge'] = topo['charge']
            atom['chargeB'] = topo['chargeB']
            atom['mass'] = topo['mass']
            atom['massB'] = topo['massB']
        else:
            # Atoms in GRO but not in topology (e.g., solvent)
            atom['type'] = 'UNK'
            atom['typeB'] = None
            atom['charge'] = 0.0
            atom['chargeB'] = None
            atom['mass'] = 0.0
            atom['massB'] = None

        merged.append(atom)

    return merged


def filter_lambda_0(atoms):
    """
    Filter atoms for lambda_0 (state A).

    Rules:
    - Include atoms where type is NOT "DUM_*"
    - Exclude atoms where type starts with "DUM_"
    """
    filtered = []

    for atom in atoms:
        atom_type = atom.get('type', 'UNK')

        # Exclude dummy atoms in state A
        if atom_type.startswith('DUM_'):
            continue

        filtered.append(atom)

    return filtered


def filter_lambda_1(atoms):
    """
    Filter atoms for lambda_1 (state B).

    Rules:
    - If typeB exists: Include atoms where typeB is NOT "DUM_*"
    - If typeB does not exist: Include atoms where type is NOT "DUM_*"
    """
    filtered = []

    for atom in atoms:
        typeB = atom.get('typeB')
        atom_type = atom.get('type', 'UNK')

        if typeB is not None:
            # Has dual topology - check state B type
            if typeB.startswith('DUM_'):
                continue  # Exclude dummy in state B
        else:
            # No dual topology - include as-is
            pass

        filtered.append(atom)

    return filtered


def write_gro(atoms, filepath, title):
    """Write atoms in GRO format"""
    with open(filepath, 'w') as f:
        f.write(f"{title}\n")
        f.write(f"{len(atoms)}\n")

        for i, atom in enumerate(atoms, 1):
            # GRO format: resnr(5) resname(5) atomname(5) atomnum(5) x(8) y(8) z(8) [vx(8) vy(8) vz(8)]
            line = f"{atom['resnr']:5d}{atom['resname']:<5s}{atom['atomname']:<5s}{atom['index']:5d}"
            line += f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"

            if atom['vx'] is not None:
                line += f"{atom['vx']:8.4f}{atom['vy']:8.4f}{atom['vz']:8.4f}"

            f.write(f"{line}\n")

        # Last line: box vectors (for now, zeros)
        f.write("0.0 0.0 0.0\n")


def build_chain_map(atoms, chain_molecules=None, topology_atoms=None, top_file=None):
    """
    Build a mapping from residue numbers to chain IDs.

    If chain_molecules provided, uses actual chain structure from topology.
    Otherwise, detects chain boundaries from residue number jumps.

    Returns dict: {resnr -> chain_id}
    """
    chain_map = {}

    if chain_molecules and top_file:
        # Use actual chain structure from topology
        from .topology_parser import extract_moleculetype_atoms
        from pathlib import Path

        chain_id_letter = 0
        for mol_name, mol_count in chain_molecules.items():
            # Find .itp file for this molecule
            top_dir = Path(top_file).parent
            itp_files = list(top_dir.glob(f"*{mol_name}*.itp"))

            if itp_files:
                residues = extract_moleculetype_atoms(str(itp_files[0]))
                chain_letter = chr(ord('A') + (chain_id_letter % 26))
                for resnr in residues:
                    chain_map[resnr] = chain_letter
                chain_id_letter += 1

    if not chain_map:
        # Fallback: detect chain boundaries from residue number jumps
        current_chain = 'A'
        prev_resnr = -1
        chain_ord = 0

        sorted_atoms = sorted(atoms, key=lambda a: a['resnr'])

        for atom in sorted_atoms:
            resnr = atom['resnr']

            if resnr < prev_resnr:
                chain_ord += 1
                if chain_ord < 26:
                    current_chain = chr(ord('A') + chain_ord)
                else:
                    current_chain = 'A'

            if resnr not in chain_map:
                chain_map[resnr] = current_chain

            prev_resnr = resnr

    return chain_map


def write_pdb(atoms, filepath, title, chain_map=None, mutation_maps=None, lambda_state=0,
              topology_atoms=None, resnr_to_mutation=None):
    """
    Write atoms in PDB format using BioPython PDBIO.

    BioPython PDBIO handles all PDB format specifications correctly including:
    - Proper column alignment
    - Atom name formatting
    - Chain IDs
    - Residue numbering

    For lambda_state=1, applies atom name renamings ONLY for atoms with dual topology
    and only uses the mutation mapping for that specific residue.
    """
    if chain_map is None:
        chain_map = build_chain_map(atoms)

    if mutation_maps is None:
        mutation_maps = {}

    if topology_atoms is None:
        topology_atoms = {}

    if resnr_to_mutation is None:
        resnr_to_mutation = {}

    # Create BioPython structure
    structure = Structure.Structure('model')
    model = Model.Model(0)
    structure.add(model)

    # Group atoms by chain and residue
    chains_dict = {}
    residues_dict = {}

    for atom in atoms:
        resnr = atom['resnr']
        resname = atom['resname']
        atom_name = atom['atomname']
        chain_id = chain_map.get(resnr, 'A')

        # For lambda_1, rename atoms that have mappings in the mutation file
        renamed_atom_name = atom_name
        if lambda_state == 1 and mutation_maps and resnr_to_mutation:
            # Get the specific mutation for this residue (e.g., 'W2Y')
            if resnr in resnr_to_mutation:
                mutation_name = resnr_to_mutation[resnr]
                if mutation_name in mutation_maps:
                    # Check if this atom has a mapping in the mutation file
                    rename_map = mutation_maps[mutation_name].get('atom_A_to_B', {})
                    if atom_name in rename_map:
                        # Apply the renaming from the mutation mapping
                        renamed_atom_name = rename_map[atom_name]
                    elif atom_name.startswith('D') and len(atom_name) > 1:
                        # Convention for dual topology atoms: D prefix indicates only in state B
                        # e.g., DHD2 is the state B version of HV -> HD2
                        # Try stripping 'D' prefix - the result should be the state B name
                        stripped_name = atom_name[1:]
                        # First check if the stripped name matches any state B atom name
                        atom_B_to_A = mutation_maps[mutation_name].get('atom_B_to_A', {})
                        if stripped_name in atom_B_to_A:
                            # This is a state B atom name, use it directly
                            renamed_atom_name = stripped_name

        # Convert from nm (GRO) to Angstroms (PDB): 1 nm = 10 Å
        coords = np.array([atom['x'] * 10.0, atom['y'] * 10.0, atom['z'] * 10.0])

        # Create/get chain
        if chain_id not in chains_dict:
            chain = Chain.Chain(chain_id)
            model.add(chain)
            chains_dict[chain_id] = chain
        else:
            chain = chains_dict[chain_id]

        # Create/get residue
        res_key = (chain_id, resnr)
        if res_key not in residues_dict:
            residue = Residue.Residue((' ', resnr, ' '), resname, 0)
            chain.add(residue)
            residues_dict[res_key] = residue
        else:
            residue = residues_dict[res_key]

        # Determine element symbol (first alphabetic character, uppercase)
        element = ''.join([c for c in renamed_atom_name if c.isalpha()])
        if element:
            element = element[0].upper()
        else:
            element = 'X'

        # Create atom and add to residue (use renamed atom name)
        # Atom(name, coord, bfactor, occupancy, altloc, fullname, serial_number, element)
        atom_obj = Atom.Atom(
            name=renamed_atom_name,
            coord=coords,
            bfactor=0.0,
            occupancy=1.0,
            altloc=' ',
            fullname=renamed_atom_name,
            serial_number=0,
            element=element
        )
        residue.add(atom_obj)

    # Write PDB file using BioPython PDBIO
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)

    # Write to temporary file first
    temp_filepath = str(filepath) + '.tmp'
    pdb_io.save(temp_filepath)

    # Read and prepend REMARK header
    with open(temp_filepath, 'r') as f:
        content = f.read()

    with open(filepath, 'w') as f:
        f.write(f"REMARK   {title}\n")
        f.write(content)

    # Remove temporary file
    import os
    os.remove(temp_filepath)


def update_ndx(lambda_0_atoms, lambda_1_atoms, input_ndx, output_ndx, residue_types=None):
    """
    Update or create index file with lambda and protein-only groups.

    Args:
        lambda_0_atoms: List of atoms in lambda_0 state
        lambda_1_atoms: List of atoms in lambda_1 state
        input_ndx: Path to input index file or None (to create new)
        output_ndx: Path to output index file
        residue_types: Dict mapping resnr to residue type ('protein', 'solvent', 'ion')
    """
    # Read existing groups (if input_ndx provided)
    groups = {}
    if input_ndx:
        try:
            with open(input_ndx) as f:
                current_group = None
                for line in f:
                    line = line.rstrip()
                    if line.startswith('['):
                        current_group = line.strip('[] ').strip()
                        groups[current_group] = []
                    elif current_group and line.strip():
                        # Parse atom indices
                        for atom_id in line.split():
                            try:
                                groups[current_group].append(int(atom_id))
                            except ValueError:
                                pass
        except FileNotFoundError:
            pass

    # Add lambda groups
    lambda_0_indices = [a['index'] for a in lambda_0_atoms]
    lambda_1_indices = [a['index'] for a in lambda_1_atoms]

    groups['lambda_0'] = lambda_0_indices
    groups['lambda_1'] = lambda_1_indices

    # Add protein-only groups (without water and ions) if residue_types provided
    if residue_types:
        lambda_0_protein = [a['index'] for a in lambda_0_atoms
                           if residue_types.get(a['resnr'], 'unknown') == 'protein']
        lambda_1_protein = [a['index'] for a in lambda_1_atoms
                           if residue_types.get(a['resnr'], 'unknown') == 'protein']

        groups['lambda_0_wo_water_and_ions'] = lambda_0_protein
        groups['lambda_1_wo_water_and_ions'] = lambda_1_protein

    # Write updated NDX file
    with open(output_ndx, 'w') as f:
        for group_name, atom_indices in groups.items():
            f.write(f"[ {group_name} ]\n")

            # Write indices 15 per line
            for i, idx in enumerate(atom_indices):
                if i > 0 and i % 15 == 0:
                    f.write("\n")
                f.write(f"{idx:6d} ")

            if atom_indices:
                f.write("\n")
