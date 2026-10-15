# -*- coding: utf-8 -*-

"""
viscosity_step
A SEAMM plug-in for Viscosity
"""

# Bring up the classes so that they appear to be directly in
# the viscosity_step package.

from .viscosity import Viscosity  # noqa: F401, E501
from .viscosity_parameters import ViscosityParameters  # noqa: F401, E501
from .viscosity_step import ViscosityStep  # noqa: F401, E501
from .tk_viscosity import TkViscosity  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
