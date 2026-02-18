"""Force field detection from topology files"""

from pathlib import Path
import re


def detect_force_field(top_file, ff_dir):
    """
    Detect force field from topology #include directives and mutres.mtp files.

    ff_dir can be either:
    - A directory containing .ff subdirectories (standard case)
    - A .ff directory itself (if user passes a .ff dir directly)

    Returns dict with 'name', 'path', 'mutres_path'
    """
    ff_info = {'name': 'unknown', 'path': None, 'mutres_path': None}

    ff_dir = Path(ff_dir)

    # If ff_dir is a .ff directory itself, use its parent
    if ff_dir.name.endswith('.ff'):
        ff_search_dir = ff_dir.parent
        ff_dir_name = ff_dir.name
    else:
        ff_search_dir = ff_dir
        ff_dir_name = None

    # Read topology file and look for #include directives
    with open(top_file) as f:
        content = f.read()

    # Look for force field includes (e.g., #include "amber99sb-star-ildn-mut.ff/forcefield.itp")
    includes = re.findall(r'#include\s*["\']([^"\']+)["\']', content)

    for include in includes:
        # Extract force field name
        if '.ff/' in include:
            # Extract the part after .ff/ to find the FF directory name
            # Handle both absolute paths and relative paths
            ff_part = include.split('.ff/')[0]
            # Get just the directory name (last component)
            ff_name = Path(ff_part).name + '.ff' if '/' in ff_part or '\\' in ff_part else ff_part + '.ff'

            # If we know the .ff dir, check if it matches
            if ff_dir_name and ff_name != ff_dir_name:
                continue

            # Search for this FF in ff_search_dir
            ff_path = ff_search_dir / ff_name
            if ff_path.exists():
                ff_info['name'] = ff_name.replace('.ff', '')
                ff_info['path'] = str(ff_path)

                # Look for mutres.mtp
                mutres_path = ff_path / 'mutres.mtp'
                if mutres_path.exists():
                    ff_info['mutres_path'] = str(mutres_path)
                break

    return ff_info
