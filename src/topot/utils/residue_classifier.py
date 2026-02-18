"""Classify residues as protein, solvent, or ion"""

from pathlib import Path


PROTEIN_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYX',  # Standard amino acids
    'GLN', 'GLU', 'GLY', 'HIS', 'HSD', 'HSE', 'HSP',  # Histidine variants
    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE',  # Selenomethionine
    'W2Y', 'W2F', 'W2L', 'W2M', 'W2A',  # Common mutations
}

SOLVENT_RESIDUES = {
    'SOL', 'WAT', 'HOH', 'H2O',  # Water
    'TIP3P', 'TIP4P', 'TIP5P', 'SPC', 'SPE',  # Water models
    'TIP3F', 'TIP3E', 'TIP4EW', 'OPC', 'OPC3',
}

ION_RESIDUES = {
    'NA', 'Na+', 'NAJ', 'NAS',  # Sodium
    'CL', 'Cl-', 'CLJ', 'CLS',  # Chloride
    'K', 'K+', 'KJ', 'KS',  # Potassium
    'MG', 'Mg2+', 'MNG',  # Magnesium
    'CA', 'Ca2+',  # Calcium
    'ZN', 'Zn2+',  # Zinc
    'LI', 'Li+',  # Lithium
    'RB', 'Rb+',  # Rubidium
}


def classify_residues(gro_atoms, ff_dir):
    """
    Classify each residue as 'protein', 'solvent', or 'ion'.

    Returns dict: {resnr -> 'protein' | 'solvent' | 'ion' | 'unknown'}
    """
    residue_types = {}

    # Try to load residuetypes.dat from FF directory
    ff_residue_types = load_residuetypes(ff_dir)

    # Classify each unique residue from GRO
    seen_residues = set()
    for atom in gro_atoms:
        resnr = atom['resnr']
        resname = atom['resname']

        if resnr in residue_types:
            continue

        # Check FF residuetypes.dat first
        if ff_residue_types and resname in ff_residue_types:
            residue_types[resnr] = ff_residue_types[resname].lower()
            continue

        # Fall back to hardcoded lists
        if resname in PROTEIN_RESIDUES:
            residue_types[resnr] = 'protein'
        elif resname in SOLVENT_RESIDUES:
            residue_types[resnr] = 'solvent'
        elif resname in ION_RESIDUES:
            residue_types[resnr] = 'ion'
        else:
            residue_types[resnr] = 'unknown'

    return residue_types


def load_residuetypes(ff_dir):
    """Load residuetypes.dat from force field directory"""
    ff_dir = Path(ff_dir)

    # Try to find residuetypes.dat
    residuetypes_file = ff_dir / 'residuetypes.dat'
    if not residuetypes_file.exists():
        return {}

    residue_types = {}
    try:
        with open(residuetypes_file) as f:
            for line in f:
                # Remove comments
                comment_idx = line.find('#')
                if comment_idx >= 0:
                    line = line[:comment_idx]

                line = line.strip()
                if not line:
                    continue

                # Split by whitespace (handles tabs and spaces)
                parts = line.split()
                if len(parts) >= 2:
                    resname = parts[0].strip()
                    restype = parts[1].strip()
                    residue_types[resname] = restype

    except IOError:
        pass

    return residue_types


def identify_mutations(topology_atoms, gro_atoms):
    """
    Identify mutated residues by comparing atoms.

    Returns list of dicts: [{
        'resnr': int,
        'atoms_disappear': int,
        'atoms_appear': int
    }, ...]
    """
    mutations = []

    # Group topology atoms by residue
    topology_by_res = {}
    for idx, atom in topology_atoms.items():
        resnr = atom['resnr']
        if resnr not in topology_by_res:
            topology_by_res[resnr] = []
        topology_by_res[resnr].append(atom)

    # Check for residues with dual topology
    for resnr, atoms in topology_by_res.items():
        # Count atoms with typeB (dual topology)
        dual_atoms = [a for a in atoms if a['typeB']]

        if not dual_atoms:
            continue  # Not a mutated residue

        # Count disappearing and appearing atoms
        disappearing = 0
        appearing = 0

        for atom in dual_atoms:
            # Atom disappears if typeB is dummy
            if atom['typeB'].startswith('DUM_'):
                disappearing += 1
            # Atom appears if type is dummy but typeB is not
            if atom['type'].startswith('DUM_') and not atom['typeB'].startswith('DUM_'):
                appearing += 1

        if disappearing > 0 or appearing > 0:
            mutations.append({
                'resnr': resnr,
                'atoms_disappear': disappearing,
                'atoms_appear': appearing
            })

    return mutations
