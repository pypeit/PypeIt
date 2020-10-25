"""
Version checking.
"""

import pkg_resources

requirements_file = pkg_resources.resource_filename('pypeit', 'requirements.txt')
install_requires = [line.strip().replace('==', '>=') for line in open(requirements_file)
                    if not line.strip().startswith('#') and line.strip() != '']
for requirement in install_requires:
    pkg, version = requirement.split('>=')
    try:
        pv = pkg_resources.get_distribution(pkg).version
    except pkg_resources.DistributionNotFound:
        raise ImportError("Package: {:s} not installed!".format(pkg))
    else:
        if pkg_resources.parse_version(pv) < pkg_resources.parse_version(version):
            raise ImportError('Your version of {0} ({1}) is incompatible with PypeIt.  '.format(pkg, pv)
                                + 'Please update to version >= {0}'.format(version))

