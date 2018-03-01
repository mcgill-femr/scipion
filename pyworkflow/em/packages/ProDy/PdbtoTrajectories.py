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

import math
import numpy as np

from pyworkflow.em import *
import pyworkflow.protocol.constants as pwconst
from prody import *
from pyworkflow.utils import *
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils.path import createLink
import commands
import matplotlib.pyplot as plt
from shutil import copyfile

class computePdbTrajectories(EMProtocol):
    """ Protocol to execute functions from ProDy software"""

    _label = "Generate Trajectories ProDy"

    def _defineParams(self, form):
        form.addSection(label="Prody Trajectories")
        form.addParam('initialPdb', PointerParam, pointerClass='PdbFile',
                      label='Initial Pdb', important=True,
                      help='Choose an initial conformation using a Pdb to '
                           'compute its N trajectories')
        form.addParam('useFinalPdb', params.BooleanParam, default=False,
                      label="Use final conformation")
        form.addParam('finalPdb', PointerParam, pointerClass='PdbFile',
                      label='Final Pdb', important=True,
                      condition = 'useFinalPdb == True',
                      help='Choose a final conformation using a Pdb to '
                           'compute its N trajectories')
        form.addParam('cycleNumber', params.IntParam, default=5,
                      label='Number of cycles',
                      condition = 'useFinalPdb == False',
                      important=True, help='Number of cycles to cover a '
                                           'landscape')
        form.addParam('numTrajectories', params.IntParam, default=5,
                      label='Number of trajectories',
                      important=True, help='Number of trajectories to generate')
        form.addParam('inputTopology', params.MultiPointerParam,
                      pointerClass='SetOfVolumes,Volume',
                      label="Input topology files", important=True,
                      help='Select topology files')
        form.addParam('inputParameters', params.MultiPointerParam,
                      label="Input parameters files", important=True,
                      help='Select parameters files')
        form.addParam('anmCutoff', params.IntParam, default=15,
                      label='ANM Cutoff', expertLevel = pwconst.LEVEL_ADVANCED,
                      help='maximum distance that two residues are in contact.')
        form.addParam('scaling', params.FloatParam, default=0.1,
                      label='Scaling (A)', expertLevel =
                      pwconst.LEVEL_ADVANCED,
                      help='The scaling factor used when disturbing the '
                           'protein structure in ANM-MC steps.')
        form.addParam('acceptanceParam', params.FloatParam, default=0.1,
                      label='Acceptance Parameter', expertLevel =
                      pwconst.LEVEL_ADVANCED,
                      help='starting value for the acceptance parameter in '
                           'ANM-MC steps')
        form.addParam('maxNumSteps', params.IntParam, default=1000000,
                      label='Max number of ANM steps', expertLevel =
                      pwconst.LEVEL_ADVANCED,
                      help='maximal number of steps in ANM-MC step.')
        form.addParam('spring', params.IntParam, default=15,
                      label='Max number of ANM steps', expertLevel=
                      pwconst.LEVEL_ADVANCED,
                      help='In targeted molecular dynamics simulation, '
                           'the target potential is harmonic and the spring '
                           'constant term shows the force applied to a given '
                           'structure to reach the target structure')


    def _insertAllSteps(self):
        self._insertFunctionStep('createTrajectories')

    def createTrajectories(self):

        topo1 = '/home/javiermota/ProDy/files_topo_param/top_all36_prot.rtf'
        topo2 = '/home/javiermota/ProDy/files_topo_param' \
                '/toppar_water_ions_mod2.str'
        param1 = '/home/javiermota/ProDy/files_topo_param/par_all36_prot.prm'
        param2 = '/home/javiermota/ProDy/files_topo_param/par_all36m_prot.prm'

        self._params = {'initPdb': self.initialPdb.get(),
                        'numTraj': self.numTrajectories.get(),
                        'cutoff': self.anmCutoff.get(),
                        'scale': self.scaling.get(),
                        'acceptParam': self.acceptanceParam.get(),
                        'maxSteps': self.maxNumSteps.get(),
                        'spr': self.spring.get()
                        }

        if self.useFinalPdb.get() is True:
            self._params['finPdb'] = self.finalPdb.get()
        else:
            self._params['cycle'] = self.cycleNumber.get()









