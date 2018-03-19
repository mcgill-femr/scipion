# **************************************************************************
# *
# * Authors:     Javier Mota (jmota@cnb.csic.es)
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
This sub-package contains data and protocol classes
wrapping Prody programs
"""
import imp

from bibtex import _bibtex # Load bibtex dict with references
from protocol_ProDy import ProdyProt
from pdbtoTrajectories import computePdbTrajectories
from computePCATrajectory import computeModesPcaPdb
from joinPDBs import joinPDBs
from viewer import ProdyViewer
from viewerPCA import ProdyViewerPca
from viewerTrajectories import ProdyTrajectoriesViewer
from protocol_import import ProtImportTrajectories
from protocol_cluster import clusterPdbTrajectories
# from viewerCluster import ProdyViewerCluster
from viewerClusterForm import ProdyViewerCluster

def validateInstallation():
    """ This function will be used to check if Prody is properly installed. """
    try:
        imp.find_module('ProDy')
        found = True
    except ImportError:
        found = False

    errors = []
    if not found:
        errors.append("ProDy not found in the system")

    return errors
