"""Parse mutres.mtp files to extract atom name mappings between states"""

from pathlib import Path


def parse_mutres_file(mutres_path):
    """
    Parse mutres.mtp file and extract atom name mappings.

    Format:
    [ mutation_name ]
    [ morphes ]
        name_A  type_A -> name_B  type_B

    Returns dict: {
        'mutation_name': {
            'atom_A_to_B': {name_A: name_B},  # State A atom name -> State B atom name
            'atom_B_to_A': {name_B: name_A}   # State B atom name -> State A atom name
        }
    }
    """
    mutations = {}

    if not Path(mutres_path).exists():
        return mutations

    # Standard GROMACS topology sections to skip
    known_sections = {
        'morphes', 'atoms', 'bonds', 'constraints', 'pairs',
        'angles', 'dihedrals', 'impropers', 'cmap',
        'coords', 'rotations'
    }

    with open(mutres_path) as f:
        lines = f.readlines()

    current_mutation = None
    in_morphes = False

    for line in lines:
        line = line.rstrip()

        # Check for section header
        if line.strip().startswith('[ ') and line.strip().find(' ]') > 0:
            # Extract section name (could have comments after ])
            section_end = line.find(']')
            section = line[line.find('[')+1:section_end].strip()

            # If it's a morphes section, mark that we're in morphes for the current mutation
            if section == 'morphes':
                in_morphes = True
                continue

            # If it's a known section but not morphes, turn off morphes parsing
            if section in known_sections:
                in_morphes = False
                continue

            # If it's not a known section, it's a mutation name
            if section not in known_sections:
                current_mutation = section
                mutations[current_mutation] = {'atom_A_to_B': {}, 'atom_B_to_A': {}}
                in_morphes = False
            continue

        # Parse morphes section
        if in_morphes and current_mutation and line.strip() and not line.strip().startswith(';'):
            # Format: name_A type_A -> name_B type_B
            parts = line.split('->')
            if len(parts) == 2:
                left = parts[0].split()
                right = parts[1].split()

                if len(left) >= 1 and len(right) >= 1:
                    name_A = left[0]
                    name_B = right[0]

                    # Skip .gone atoms - these are virtual markers, not real atoms
                    if '.gone' in name_A or '.gone' in name_B:
                        continue

                    # Only store mapping if both atoms could exist (no DUM_ prefix in base names)
                    # Store mappings both ways
                    mutations[current_mutation]['atom_A_to_B'][name_A] = name_B
                    mutations[current_mutation]['atom_B_to_A'][name_B] = name_A

    return mutations


def get_atom_rename_map(mutation_name, mutres_path, direction='B_to_A'):
    """
    Get atom name mapping for a specific mutation.

    Args:
        mutation_name: Name of mutation (e.g., 'W2Y')
        mutres_path: Path to mutres.mtp file
        direction: 'A_to_B' (state A names to state B) or 'B_to_A' (state B names to state A)

    Returns:
        dict: {old_name -> new_name}
    """
    mutations = parse_mutres_file(mutres_path)

    if mutation_name not in mutations:
        return {}

    if direction == 'A_to_B':
        return mutations[mutation_name]['atom_A_to_B']
    else:  # B_to_A
        return mutations[mutation_name]['atom_B_to_A']
