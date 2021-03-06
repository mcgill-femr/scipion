# **************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains the protocols and data for LocScale
"""
from locscale import *
from bibtex import _bibtex # Load bibtex dict with references

_logo = "locscale_logo.jpg"

from protocol_locscale import ProtLocScale

def validateInstallation():
    """ This function will be used to check if package is properly installed."""
    missingPaths = ["%s: %s" % (var, os.environ[var])
                    for var in [LOCSCALE_HOME_VAR, EMAN2DIR_VAR]
                    if not os.path.exists(os.environ[var])]

    if missingPaths:
        return ["Required software not found in the system:"] + missingPaths + \
               ["Try to install the package following: scipion install --help"]
    else:
        return [] # No errors