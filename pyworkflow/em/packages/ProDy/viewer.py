# **************************************************************************
# *
# * Authors:     Javier Mota Garcia (jmota@cnb.csic.es)
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
This module implement the wrappers around ProDy protocol
visualization program.
"""

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from pyworkflow.gui.text import *
from pyworkflow.gui.dialog import showError, showWarning
from pyworkflow.gui.plotter import Plotter
import glob
from protocol_ProDy import ProdyProt
from prody import *


class ProdyViewer(Viewer):
    """ Wrapper to visualize Pdb to SAXS. """
    _targets = [ProdyProt]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, **args):
        Viewer.__init__(self, **args)

    def visualize(self, obj, **args):

        cls = type(obj)

        if issubclass(cls, ProdyProt):
            paths = open(obj._getExtraPath("paths.txt"), "r")
            content = paths.readlines()
            content = [line.rstrip('\n') for line in content]

            fnAnm = obj._getExtraPath("anmModes.anm.npz")
            fnPca = obj._getExtraPath("pcaModes.pca.npz")
            initPdb = content[0]
            fnInitTraj = content[1]
            fnFinTraj = content[2]

            Pdb = parsePDB(initPdb)

            anm = loadModel(fnAnm)
            pca = loadModel(fnPca)

            showOverlapTable(pca, anm)

            initTraj = Trajectory(fnInitTraj)
            initTraj.setCoords(Pdb) # Set the initial structure as the reference
            initTraj.setAtoms(Pdb.ca)  # A shortcut for .select('ca')

            finTraj = Trajectory(fnFinTraj)
            finTraj.setCoords(Pdb)  # Set the initial structure as the reference
            finTraj.setAtoms(Pdb.ca)  # A shortcut for .select('ca')

            showProjection(initTraj, pca[:2], color='g', new_fig=True)
            showProjection(finTraj, pca[:2], color='r',  new_fig=False)

            showProjection(initTraj, anm[:3], color='g', new_fig=True)
            showProjection(finTraj, anm[:3], color='r',  new_fig=False)



