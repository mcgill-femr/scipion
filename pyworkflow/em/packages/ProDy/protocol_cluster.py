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

class clusterPdbTrajectories(EMProtocol):
    """ Protocol to execute functions from ProDy software"""

    _label = "Cluster Trajectories ProDy"

    def _defineParams(self, form):
        self.defaultCycles = 5

        form.addSection(label="ProDy Cluster Trajectories")
        form.addParam('setOfTrajectories', PointerParam,
                      pointerClass='SetOfTrajectories',
                      label="Set of Trajectories",
                      important=True)
        form.addParam('distanceMatrixFile', PointerParam,
                      pointerClass='EMFile',
                      label="Distance Matrix File",
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

        ens = parseDCD(str(fileName))
        pdb = parsePDB(str(trajectory._initialPdb))
        ca = pdb.ca.copy()
        ens.setCoords(ca)
        ens.setAtoms(ca)

        fnPdb = []
        setOfPDBs = self._createSetOfPDBs()
        numbers = []
        for i, conformation in enumerate(ens):
            fnPdb.append(self._getExtraPath('pdb{:02d}.pdb'.format(i)))
            writePDB(fnPdb[i], conformation)
            pdb = PdbFile(fnPdb[i])
            setOfPDBs.append(pdb)
            numbers.append(str(i))

        distanceMatrix = parseArray(str(self.distanceMatrixFile
                                        .get().getFileName()))

        self.tree = calcTree(names=numbers, distance_matrix=distanceMatrix,
                             method=treeMethod)
        self.subgroups = findSubgroups(self.tree, self.subgroupCutoff.get())

        self.setOfRepresentatives = self._createSetOfPDBs()
        for subgroup in self.subgroups:
            print subgroup

            repName = subgroup[0]
            for i in subgroup:
                distList = i, mean([distanceMatrix[int(i)][int(j)]
                                    for j in subgroup])
                if distList[-1] < minDist:
                    minDist = distList[-1]
                    repName = distList[0]

            print repName, minDist

            repPdb = PdbFile(self._getExtraPath('pdb{:02d}.pdb'
                                             .format(int(repName))))

            self.setOfRepresentatives.append(repPdb)

    def createOutputStep(self):
        Phylo.write(self.tree, self._getExtraPath('clustering_tree.nwk'),
        'newick')
        myFile = EMFile(self._getExtraPath('clustering_tree.nwk'))
        self._defineOutputs(ClusteringTree=myFile)

        self._defineOutputs(representivePDBs=self.setOfRepresentatives)