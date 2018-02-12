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
from prody import *
from pyworkflow.em.packages.ProDy import PRODY
from pyworkflow.utils import *
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils.path import createLink
import commands
import matplotlib.pyplot as plt

class ProdyProt(EMProtocol):

    """ Protocol to execute functions from ProDy software"""

    _label = "ProDy protocol"

    def _defineParams(self, form):

        form.addSection(label = "Prody NMA analysis and molecular dynamics")
        form.addParam('inputStructure', PathParam, label="Input structure",
                      important=True,
                      #pointerClass='PdbFile',
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model(an EM volume '
                           'converted into pseudoatoms)')
        form.addParam('initTrajectory', PathParam,
                      label="Initial Trajectory",
                      important=True)
        form.addParam('finTrajectory', PathParam, label="Final Trajectory",
                      important=True)

        form.addSection(label='Animation')
        form.addParam('amplitude', FloatParam, default=50,
                      label="Amplitude")
        form.addParam('nframes', IntParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of frames')
        form.addParam('downsample', FloatParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Downsample pseudoatomic structure',
                      help='Downsample factor 2 means removing one half of '
                           'the atoms or pseudoatoms.')
        form.addParam('pseudoAtomThreshold', FloatParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Pseudoatom mass threshold',
                      help='Remove pseudoatoms whose mass is below this '
                           'threshold. This value should be between 0 and 1.\n'
                           'A threshold of 0 implies no atom removal.')

    def _insertAllSteps(self):
        self._insertFunctionStep('prodyWrapper')

    def prodyWrapper(self):
        time.sleep(10)
        inputStructure = self.inputStructure.get().getFileName()
        #os.path.abspath(self.inputStructure.get().getFileName())
        initialTrajectory = self.initTrajectory.get().getFileName()
        finaltrajectory = self.finTrajectory.get().getFileName()
        #createLink(inputStructure, 'pseudoatoms.pdb')
        #self.runJob(PRODY,"pseudoatoms.pdb")
        GluA2_sim_ca = parsePDB(inputStructure, subset='ca')


    def _validate(self):
        errors = []
        if which('prody') is '':
            errors.append('You should have the program crysol in the PATH')
        return errors


