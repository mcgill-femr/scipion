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
import pyworkflow.protocol.constants as pwconst

class computeModesPcaPdb(EMProtocol):

    _label = "Compute PCA Pdb trajectories"

    def _defineParams(self, form):
        form.addSection(label="Prody PCA Analysis")
        form.addParam('setType', params.EnumParam,
                      choices=['setOfTrajectories', 'setOfPDBs'],
                      default=0,
                      label="Input structures format")
        form.addParam('setOfTrajectories', PointerParam,
                      pointerClass='SetOfTrajectories',
                      label="Set of Trajectories",
                      condition='setType == 0',
                      important=True)
        form.addParam('setOfPDBs', PointerParam,
                      pointerClass='SetOfPDBs',
                      label="Set of PDBs",
                      condition='setType == 1',
                      important=True)
        form.addParam('usePdb', params.BooleanParam, default=False,
                      label="Use a separate initial PDB",
                      help = "If you select this option you must provide a "
                             "PdbFile that will be used as initial PDB for "
                             "the input set",
                      expertLevel=pwconst.LEVEL_ADVANCED)
        form.addParam('Pdb', PointerParam, pointerClass='PdbFile',
                      label='Input Pdb', condition='usePdb == True',
                      help="Select the PdbFile that will be used as initial "
                           "PDB for the input set",
                      expertLevel=pwconst.LEVEL_ADVANCED)

    def _insertAllSteps(self):
        self._insertFunctionStep('_calcPCA')
        self._insertFunctionStep('_createOutputStep')

    def _calcPCA(self):

        if self.Pdb.get() is not None:
            self.pdbFileName = self.Pdb.get().getFileName()
        else:
            foundPdb = False
            if self.setType == 0:
                for traj in self.setOfTrajectories.get():
                    if traj._initialPdb.get() is not None:
                        foundPdb = True
                        self.pdbFileName = traj._initialPdb.get()
                        break
            else:
                #pdbFileNames = []
                for i, pdb in enumerate(self.setOfPDBs.get()):
                    if i == 0:
                        foundPdb = True
                        self.pdbFileName = str(pdb.getFileName())
                    #pdbFileNames.append(str(pdb.getFileName()))

            if not foundPdb:
                self.pdbFileName = None
                raise ValueError("There needs to be at least one PDB associated"
                                 " with the setOfTrajectories or setOfPDBs or "
                                 "provided separately.")

        if self.pdbFileName.endswith('.pdb'):
            pdb = parsePDB(self.pdbFileName)
        elif self.pdbFileName.endswith('.cif'):
            pdb = parseCIF(self.pdbFileName)
        else:
            raise ValueError('The starting PDB filename needs to end in '
                             '.pdb or .cif')

        if self.setType == 1:

            #TODO: adapt this part when the input is a setofpdbs to pseudoatoms
            pdbs = []
            for pdbFn in self.setOfPDBs.get():
                if pdbFn.endswith('.pdb'):
                    pdbs.append(parsePDB(pdbFn,subset='ca'))
                elif pdbFn.endswith('.cif'):
                    pdbs.append(parseCIF(pdbFn,subset='ca'))
                else:
                    raise ValueError('All PDB filenames should end in'
                                     ' .pdb or .cif')

            for i, currPdb in enumerate(pdbs):
                #flag = currPdb.getPseudoAtoms()
                if currPdb.numAtoms() != pdbs[0].numAtoms():
                    pdbs.pop(i)
                    print("Scipion can only use sets of PDBs where all "
                          "members have the same number of atoms. Therefore"
                          " %s was removed as it does not have the same "
                          "number of atoms as the reference PDB."
                          %currPdb.getTitle())

            self.ens = buildPDBEnsemble(pdb,pdbs,mapping_func=mapChainByChain)

        else:
            combined_traj = Trajectory("combined traj for PCA")
            for traj in self.setOfTrajectories.get():
                self.flag = traj.isPseudoatoms()
                combined_traj.addFile(traj.getFileName())
            if self.flag:
                combined_traj.setCoords(pdb)
                combined_traj.setAtoms(pdb)

                self.ens = Ensemble(combined_traj)
                self.ens.setCoords(pdb)

            else:
                combined_traj.setCoords(pdb)
                combined_traj.setAtoms(pdb.ca)

                self.ens = Ensemble(combined_traj)
                self.ens.setCoords(pdb.ca)
            combined_traj.getCoordsets()
            for i, coordset in enumerate(combined_traj.getCoordsets()):
                self.ens.addCoordset(coordset)

        if self.flag:
            self.ens.setAtoms(pdb)
        else:
            self.ens.setAtoms(pdb.ca)

        self.ens.superpose()

        self.pca = PCA('all trajectories')
        self.pca.buildCovariance(self.ens)
        self.pca.calcModes()

        projection = calcProjection(self.ens, self.pca[:2], norm=False)

        self.distanceMatrix = zeros((len(self.ens),len(self.ens)))

        for i in range(len(self.ens)):
            for j in range(len(self.ens)):
                for k in range(2):
                    self.distanceMatrix[i][j] += \
                    (projection[i][k] - projection[j][k])**2
                self.distanceMatrix[i][j] = sqrt(self.distanceMatrix[i][j])


    def _createOutputStep(self):

        if self.setType.get() == 1:
            sigma = self.setOfPDBs.get().getDeviation()
        else:
            sigma = self.setOfTrajectories.get().getDeviation()

        writeDCD(self._getExtraPath('combined_trajectory.dcd'),self.ens)
        setOfTrajectories = self._createSetOfTrajectories()
        trajectory = TrajectoryDcd(self._getExtraPath('combined_trajectory.dcd')
                                   , self.pdbFileName, self.flag, sigma=sigma)
        setOfTrajectories.append(trajectory)
        self._defineOutputs(combinedTrajectory=setOfTrajectories)

        fnOutPca = self._getExtraPath("pcaModes")
        saveModel(self.pca, fnOutPca)
        myFile = EMFile(fnOutPca)
        self._defineOutputs(pcaNpzFile=myFile)

        writeArray(self._getExtraPath('distance_matrix.txt'),
                   self.distanceMatrix,'%.18e','\t')
        myFile = EMFile(self._getExtraPath('distance_matrix.txt'))
        self._defineOutputs(distanceMatrixFile=myFile)

