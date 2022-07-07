"""The testsystems used by Transformato"""

# Add imports here
from .transformato_testsystems import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
