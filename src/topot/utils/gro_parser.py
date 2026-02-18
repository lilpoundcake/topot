"""Parse GROMACS GRO coordinate files"""


def parse_gro(gro_file):
    """
    Parse GRO file and extract atomic coordinates.

    Returns list of dicts: [{
        'index': int,        # Atom number from GRO
        'resnr': int,
        'resname': str,
        'atomname': str,
        'x': float,
        'y': float,
        'z': float,
        'vx': float | None,
        'vy': float | None,
        'vz': float | None
    }, ...]
    """
    atoms = []

    with open(gro_file) as f:
        lines = f.readlines()

    # Line 1: title
    # Line 2: number of atoms
    if len(lines) < 2:
        return atoms

    try:
        n_atoms = int(lines[1].strip())
    except ValueError:
        return atoms

    # Parse atom lines (lines 3 to 3+n_atoms)
    for i in range(2, min(2 + n_atoms, len(lines))):
        line = lines[i]

        # GRO format (fixed width):
        # residue number (5 positions, integer)
        # residue name (5 positions, character string)
        # atom name (5 positions, character string)
        # atom number (5 positions, integer)
        # x, y, z, vx, vy, vz (8 positions each, real)

        try:
            resnr = int(line[0:5].strip())
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            atom_num = int(line[15:20].strip())

            x = float(line[20:28].strip())
            y = float(line[28:36].strip())
            z = float(line[36:44].strip())

            # Velocities are optional
            vx = None
            vy = None
            vz = None
            if len(line) >= 68:
                try:
                    vx = float(line[44:52].strip())
                    vy = float(line[52:60].strip())
                    vz = float(line[60:68].strip())
                except ValueError:
                    pass

            atoms.append({
                'index': atom_num,
                'resnr': resnr,
                'resname': resname,
                'atomname': atomname,
                'x': x,
                'y': y,
                'z': z,
                'vx': vx,
                'vy': vy,
                'vz': vz
            })

        except (ValueError, IndexError):
            continue

    return atoms
