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
from prody import *
from numpy import zeros, sqrt
import time

FILE = 0
PDB = 1

class computeModesPcaPdb(EMProtocol):

    _label = "Compute PCA Pdb trajectories"

    def _defineParams(self, form):
        form.addSection(label="Prody PCA Analysis")
        form.addParam('usePdb', params.BooleanParam, default=False,
                      label="Use a separate PDB")
        form.addParam('FilePdb', params.EnumParam, choices=['File', 'Pdb'],
                      default=0, important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      condition='usePdb == True',
                      label="Input file or pdb")
        form.addParam('inputStructure', PathParam, label="Input structure",
                      important=True,
                      condition='FilePdb == %s and usePdb == True' % FILE,
                      help='The input structure can be an atomic model '
                           '(true PDB) or a pseudoatomic model (an EM volume '
                           'converted into pseudoatoms)')
        form.addParam('Pdb', PointerParam, pointerClass='PdbFile',
                      label='Input Pdb', important=True,
                      condition='FilePdb == %s and usePdb == True' % PDB)
        form.addParam('setOfTrajectories', PointerParam,
                      pointerClass='SetOfTrajectories',
                      label="Set of Trajectories",
                      important=True)


    def _insertAllSteps(self):
        self._insertFunctionStep('_calcPCA')
        self._insertFunctionStep('_createOutputStep')

    def _calcPCA(self):
        time.sleep(10)

        if self.inputStructure.get() is not None:
            self.pdbFileName = self.inputStructure.get()

        elif self.Pdb.get() is not None:
            self.pdbFileName = self.Pdb.get().getFileName()

        else:
            foundPdb = False
            for traj in self.setOfTrajectories.get():
                if traj._initialPdb.get() is not None:
                    foundPdb = True
                    self.pdbFileName = traj._initialPdb.get()
                    break

            if not foundPdb:
                self.pdbFileName = None
                raise ValueError("There needs to be at least one PDB associated "
                                 "with a trajectory or provided separately.")


        pdb = parsePDB(self.pdbFileName)
        combined_traj = Trajectory("combined traj for PCA")
        for traj in self.setOfTrajectories.get():
            combined_traj.addFile(traj.getFileName())

        combined_traj.setCoords(pdb)

        self.ens = Ensemble(combined_traj)
        self.ens.setCoords(pdb)

        for i, coordset in enumerate(combined_traj.getCoordsets()):
            self.ens.addCoordset(coordset)

        self.ens.setAtoms(pdb.ca)
        self.ens.superpose()

        self.pca = PCA('all trajectories')
        self.pca.buildCovariance(self.ens)
        self.pca.calcModes()

    def _createOutputStep(self):
        n = 2
        projection = calcProjection(self.ens, self.pca[:n],
                                    norm=False)
        self.distanceMatrix = zeros((len(self.ens),len(self.ens)))

        for i in range(len(self.ens)):
            for j in range(len(self.ens)):
                for k in range(n):
                    self.distanceMatrix[i][j] += \
                    (projection[i][k] - projection[j][k])**2

                self.distanceMatrix[i][j] = sqrt(self.distanceMatrix[i][j])

        writeArray(self._getExtraPath('distance_matrix.txt'),
                   self.distanceMatrix,'%.18e','\t')
        myFile = EMFile(self._getExtraPath('distance_matrix.txt'))
        self._defineOutputs(distanceMatrix=myFile)

        writeDCD(self._getExtraPath('combined_trajectory.dcd'),self.ens)
        setOfTrajectories = self._createSetOfTrajectories()
        trajectory = TrajectoryDcd(self._getExtraPath('combined_trajectory.dcd')
                                   , self.pdbFileName)
        setOfTrajectories.append(trajectory)
        self._defineOutputs(combinedTrajectory=setOfTrajectories)

        fnOutPca = self._getExtraPath("pcaModes")
        saveModel(self.pca, fnOutPca)
        myFile = EMFile(fnOutPca)
        self._defineOutputs(pcaNpzFile=myFile)

