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
from Bio import Phylo
from numpy import mean
from numpy import zeros, sqrt

class clusterPdbTrajectories(EMProtocol):
    """ Protocol to execute functions from ProDy software"""

    _label = "Cluster Trajectories ProDy"

    def _defineParams(self, form):

        form.addSection(label="ProDy Cluster Trajectories")
        form.addParam('setOfTrajectories', PointerParam,
                      pointerClass='SetOfTrajectories',
                      label="Set of Trajectories",
                      important=True)
        form.addParam('subgroupCutoff', params.FloatParam, default=0.01,
                      label='Distance cutoff for defining clusters',
                      expertLevel=pwconst.LEVEL_ADVANCED,
                      help='The distance in PCA space that is used for '
                           'dividing the tree into subgroups for hierarchical '
                           'clustering.')
        form.addParam('treeMethod', params.EnumParam, choices=['nj', 'upgma'],
                      expertLevel=pwconst.LEVEL_ADVANCED,
                      default=0,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label="Method to use for tree construction")
        form.addParam('pcaNpzFile', PointerParam,
                      pointerClass='EMFile',
                      label="PCA NPZ File",
                      important=True)
        form.addParam('distanceType', params.EnumParam,
                      choices=['PCA', 'RMSD'],
                      default=0,
                      label="Type of distance measure used for clustering")
        form.addParam('numModes', params.IntParam, default=2,
                      label='Number of modes',
                      help='Number of modes for visualizing '
                           '(2 or 3). This number would also '
                           'be used for distance calculations '
                           'if PCA is selected as distance '
                           'measure.',
                      condition='distanceType == 0')
        form.addParam('sampling', params.FloatParam, default=1.0,
                      label="Sampling rate (A/px)",
                      help='Sampling rate (Angstroms/pixel)')
        form.addParam('size', params.IntParam,
                      allowsNull=True,
                      label="Final size (px)",
                      help='Final size in pixels.')

    def _insertAllSteps(self):
        self._insertFunctionStep('clusterTrajectories')
        self._insertFunctionStep('createOutputStep')

    def clusterTrajectories(self):

        if self.treeMethod.get() == 0:
            treeMethod = 'nj'
        else:
            treeMethod = 'upgma'

        #import time
        #time.sleep(10)

        for trajectory in self.setOfTrajectories.get():
            fileName = trajectory.getFileName()
            break

        self.ens = parseDCD(str(fileName))

        if str(trajectory._initialPdb).endswith('.pdb'):
            pdb = parsePDB(str(trajectory._initialPdb))
        elif str(trajectory._initialPdb).endswith('.cif'):
            pdb = parseCIF(str(trajectory._initialPdb))
        else:
            raise ValueError('The initial PDB should have a filename '
                             'ending in .pdb or .cif')

        ca = pdb.ca.copy()
        self.ens.setCoords(ca)
        self.ens.setAtoms(ca)

        fnPdb = []
        setOfPDBs = self._createSetOfPDBs()
        numbers = []
        for i, conformation in enumerate(self.ens):
            fnPdb.append(self._getExtraPath('pdb{:02d}.pdb'.format(i)))
            writePDB(fnPdb[i], conformation)
            pdb = PdbFile(fnPdb[i])
            setOfPDBs.append(pdb)
            numbers.append(str(i))

        self.pca = loadModel(str(self.pcaNpzFile.get().getFileName())
                             + '.pca.npz')

        projection = calcProjection(self.ens, self.pca[:int(
            self.numModes.get())], norm=False)

        self.distanceMatrix = zeros((len(self.ens),len(self.ens)))

        for i in range(len(self.ens)):
            for j in range(len(self.ens)):
                if self.distanceType.get() == 0:
                    for k in range(int(self.numModes.get())):
                        self.distanceMatrix[i][j] += \
                        (projection[i][k] - projection[j][k])**2

                    self.distanceMatrix[i][j] = sqrt(self.distanceMatrix[i][j])
                else:
                    self.distanceMatrix[i][j] = calcRMSD(self.ens[i],
                                                         self.ens[j])

        self.tree = calcTree(names=numbers, distance_matrix=self.distanceMatrix,
                             method=treeMethod)
        self.subgroups = findSubgroups(self.tree, self.subgroupCutoff.get())

        self.setOfRepresentatives = self._createSetOfPDBs()

        print(len(self.distanceMatrix)-1)

        for subgroup in self.subgroups:
            print subgroup

            if '0' in subgroup:
                repName = '0'
            elif str(len(self.distanceMatrix) - 1) in subgroup:
                repName = str(len(self.distanceMatrix) - 1)
            else:
                minDist = self.subgroupCutoff.get()
                repName = subgroup[0]
                for i in subgroup:
                    distList = i, mean([self.distanceMatrix[int(i)][int(j)]
                                        for j in subgroup])
                    if distList[-1] < minDist:
                        minDist = distList[-1]
                        repName = distList[0]

                print repName, minDist
            print repName

            repPdb = PdbFile(self._getExtraPath('pdb{:02d}.pdb'
                                             .format(int(repName))))

            self.setOfRepresentatives.append(repPdb)

    def createOutputStep(self):

        writeArray(self._getExtraPath('distance_matrix.txt'),
                   self.distanceMatrix,'%.18e','\t')
        myFile = EMFile(self._getExtraPath('distance_matrix.txt'))
        self._defineOutputs(distanceMatrixFile=myFile)

        Phylo.write(self.tree, self._getExtraPath('clustering_tree.nwk'),
        'newick')
        myFile = EMFile(self._getExtraPath('clustering_tree.nwk'))
        self._defineOutputs(ClusteringTree=myFile)

        self._defineOutputs(representativePDBs=self.setOfRepresentatives)

        outputVols = self._createSetOfVolumes()
        outputVols.setSamplingRate(self.sampling.get())
        for i, pdb in enumerate(self.setOfRepresentatives):
            outFile = self._getExtraPath('output_vol%d'%(i+1))
            args = '-i %s --sampling %f -o %s' % (pdb.getFileName(),
                                                  self.sampling.get(), outFile)
            args += ' --size %d  --centerPDB' % self.size.get()
            program = "xmipp_volume_from_pdb"
            self.runJob(program, args)
            outVol = Volume()
            outVol.setSamplingRate(self.sampling.get())
            outVol.setFileName(outFile+'.vol')
            outputVols.append(outVol)
        outputVols.setDim(ImageDim(self.size.get(), self.size.get(), self.size.get()))

        self._defineOutputs(outputVolumes=outputVols)
