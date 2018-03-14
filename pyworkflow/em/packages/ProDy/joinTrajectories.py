# **************************************************************************
# *
# * Authors:  Javier Mota Garcia (jmota@cnb.csic.es), February 2018
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

from pyworkflow.em import *
import pyworkflow.protocol.constants as pwconst
from prody import *
from pyworkflow.utils import *
import time
import glob

class joinTrajectoriesFiles(EMProtocol):


    def _defineParams(self, form):

        form.addSection("Trajectories")
        form.addParam('trajFile', params.FileParam, label='Files directory',
                      help='Choose the trajectory files directory')
        form.addParam('pattern', params.StringParam, label = 'Pattern',
                      help = 'The pattern must be end by .dcd. If there are '
                             'more .dcd files in your directory try to put'
                             'yourword*.dcd')

    def _insertAllSteps(self):
        self._insertFunctionStep('joinFiles')

    def joinFiles(self):
        time.sleep(15)
        filesPath = self.trajFile.get('').strip()
        filesPattern = self.pattern.get('').strip()

        if filesPattern:
            fullPattern = join(filesPath, filesPattern)
        else:
            fullPattern = filesPath

        pattern = expandPattern(fullPattern.replace("$", ""))
        files = glob.glob(pattern)
        setOfTrajectories = self._createSetOfTrajectories()
        all_trajectories = Trajectory("all trajectories combined")
        for f in files:
            all_trajectories.addFile(f)
            setOfTrajectories.append(f)
        writeDCD(self._getExtraPath("trajectories_combined.dcd"),
                 all_trajectories)
        self._defineOutputs(outputTrajs=setOfTrajectories)
        self._defineSourceRelation(self.initialPdb.get(), setOfTrajectories)

