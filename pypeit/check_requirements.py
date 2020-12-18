"""
Version checking.
"""

import importlib
import pkg_resources
from pypeit import msgs

requirements_file = pkg_resources.resource_filename('pypeit', 'requirements_imports.txt')
install_requires = [line.strip().replace('==', '>=') for line in open(requirements_file)
                    if not line.strip().startswith('#') and line.strip() != '']
for requirement in install_requires:
    pkg, version = requirement.split('>=', maxsplit=1)
    version, *not_found_msg = version.split(",", maxsplit=1)
    try:
        pkg = importlib.import_module(pkg)
        pv = pkg.__version__
    except ModuleNotFoundError:
        if not_found_msg:
            msgs.warn(f"Package {pkg} not found. {not_found_msg[0]}")
        else:
            msgs.warn(f"Package {pkg} not found.")
    else:
        if pkg_resources.parse_version(pv) < pkg_resources.parse_version(version):
            if not_found_msg:
                msgs.warn(f"Your {pkg.__name__} (version {pv}) is incompatible"
                    f" with PypeIt. Please update to version >= {version}. "
                    f"{not_found_msg[0]}")
            else:
                msgs.warn(f"Your {pkg.__name__} (version {pv}) is incompatible"
                    f" with PypeIt. Please update to version >= {version}")
