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

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from pyworkflow.gui.project import ProjectWindow
from pdbtoTrajectories import computePdbTrajectories
from protocol_import import ProtImportTrajectories
from prody import *

class ProdyTrajectoriesViewer(Viewer):

    _targets = [computePdbTrajectories, ProtImportTrajectories]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, **args):
        Viewer.__init__(self, **args)

    def visualize(self, obj, **args):
        showVmdView(self.protocol)

def createVmdView(protocol):
    mystring = "-e viewTrajectories.vmd"
    fhCmd = open("viewTrajectories.vmd", 'w')

    trajectories = protocol.outputTrajs
    longestLength = 0
    longestTrajID = len(trajectories)
    for i, trajectory in enumerate(trajectories):
        fname = str(trajectory._filename)

        initPdb = trajectory.getInitialPdb()
        if initPdb == None:
            prefix = '.'.join(fname.split('.')[:-1])
            fhCmd.write("mol new " + prefix + "_pdb01.pdb\n")
        else:
            fhCmd.write("mol new " + str(initPdb) + "\n")

        fhCmd.write("mol addfile " + fname + "\n")
        fhCmd.write("animate delete  beg 0 end 0 skip 0 %d\n" %i)

        if len(Trajectory(fname)) > longestLength:
            longestLength = len(Trajectory(fname))
            longestTrajID = i

    fhCmd.write("mol top {0}\n".format(longestTrajID))
    fhCmd.close()

    return VmdView("%s" % mystring)

def showVmdView(protocol):
    createVmdView(protocol).show()
