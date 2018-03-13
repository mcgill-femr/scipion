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
from protocol_ProDy import ProdyProt
from pyworkflow.utils import *
import time
from shutil import copy

FILE = 0
PDB = 1

class computeModesPcaPdb(EMProtocol):

    _label = "Compute PCA Pdb trajectories"

    def _defineParams(self, form):
        form.addSection(label="Prody PCA Analysis")
        form.addParam('FilePdb', params.EnumParam, choices=['File', 'Pdb'],
                      default=0, important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label="Input file or pdb")
        form.addParam('inputStructure', PathParam, label="Input structure",
                      important=True,
                      condition='FilePdb == %s' % FILE,
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model (an EM volume '
                           'converted into pseudoatoms)')
        form.addParam('Pdb', PointerParam, pointerClass='PdbFile',
                      label='Input Pdb', important=True,
                      condition='FilePdb == %s' % PDB)

        form.addParam('initTrajectory', PathParam,
                      label="Initial Trajectory",
                      important=True)

        #form.addParam('finTrajectory', PathParam, label="Final Trajectory",
                      #important=True)


    def _insertAllSteps(self):
        self._insertFunctionStep('_calcPCA')

    def _calcPCA(self):
        createLink(self.inputStructure.get(), self._getExtraPath(
            'inputPdb.pdb'))
        createLink(self.initTrajectory.get(), self._getExtraPath(
            'initTraj.dcd'))
        print type(self.inputStructure.get())
        time.sleep(10)
        sim = parsePDB(self.inputStructure.get())
        protein = sim.select('name CA').copy()
        #initialTrajectory = self.initTrajectory.get()
        #finalTrajectory = self.finTrajectory.get()
        combined_traj = Trajectory(self.initTrajectory.get())
        combined_traj.setCoords(protein)
        combined_traj.setAtoms(protein)
        '''combined_traj.addFile(finalTrajectory)'''
        self.pca = PCA('AMPAR trajectories')
        self.pca.buildCovariance(combined_traj)
        self.pca.calcModes()
        fnOutPca = self._getExtraPath("pcaModes")
        saveModel(self.pca, fnOutPca)

