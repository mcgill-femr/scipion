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
from protocol_cluster import clusterPdbTrajectories
from prody import *
import matplotlib.pyplot as plt
from Bio import Phylo
from numpy import array, where
from matplotlib.pyplot import gca
from pyworkflow.protocol.params import LabelParam, BooleanParam

PDB_ALL = 0
PDB_SEL = 1
PDB_CHOICES = ['all', 'selection']


class ProdyViewerCluster(ProtocolViewer):
    """ Wrapper to visualize Pdb to SAXS. """
    _label = 'viewer ProDy clustering'
    _targets = [clusterPdbTrajectories]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, **args):
        ProtocolViewer.__init__(self, **args)


    def _defineParams(self, form):
        form.addSection(label='Visualization')
        # Select what to show
        form.addParam('doShowPdb', EnumParam, default=PDB_ALL,
                      choices=PDB_CHOICES,
                      display=EnumParam.DISPLAY_HLIST,
                      label="PDBs to visualize")
        form.addParam('showSeveralPdbs', StringParam, default='',
                      label='PDBs selection',
                      condition='doShowPdb==%d' % PDB_SEL,
                      help='Specify a list of PDBs like: 0, 1, 3 or 0-3 ')
        form.addParam('showPdb', LabelParam,
                      label="Visualize PDBs",
                      help="It shows the selected PDBs.")
        form.addParam('showTree', LabelParam,
                      label="Visualize clustering tree",
                      help="It shows the clustering tree.")
        form.addParam('showDistMatrix', LabelParam,
                      label="Visualize distance matrix",
                      help="It shows the distance matriz between PDBs.")
        form.addParam('showClusterProj', LabelParam,
                      label="Visualize projected clusters",
                      help="Visualize projected PDBs in the PCA space "
                            "representing the different clusters with "
                            "different colors")



    def _getVisualizeDict(self):
        return {'showPdb': self._viewPdb,
                'showTree': self._viewTree,
                'showDistMatrix': self._viewDistMatrix,
                'showClusterProj': self._viewClusterProj,
                }

    def _viewPdb(self, e=None):

        obj = self.protocol
        mystring = ""

        if self.doShowPdb.get() == PDB_ALL:
            for pdb in obj.representativePDBs:
                mystring += "%s " %pdb.getFileName()

        if self.doShowPdb.get() == PDB_SEL:
            nums = self.showSeveralPdbs.get()
            listNum = nums.split('-')
            ini = listNum[0]
            if len(listNum)>1:
                end = listNum[1]
            else:
                end = ini
            for i, pdb in enumerate(obj.representativePDBs):
                if i>=int(ini) and i<=int(end):
                    mystring += "%s " % pdb.getFileName()

        newVmd = VmdView("%s" % mystring)
        newVmd.show()


    def _viewTree(self, e=None):

        obj = self.protocol
        cls = type(obj)

        if issubclass(cls, clusterPdbTrajectories):

            tree = Phylo.read(str(obj.ClusteringTree.getFileName()),
                              'newick')

            subgroup_color_dict = {}
            labels = []
            for group in findSubgroups(tree, float(obj.subgroupCutoff.get())):
                c = (random.random(), random.random(), random.random())
                for label in group:
                    subgroup_color_dict[label] = c
                    labels.append(label)

            # plt.figure()
            show = showTree(tree, format='plt',
                            label_colors=subgroup_color_dict)

            plt.show()

            return show


    def _viewDistMatrix(self, e=None):

        obj = self.protocol
        cls = type(obj)

        if issubclass(cls, clusterPdbTrajectories):

            tree = Phylo.read(str(obj.ClusteringTree.getFileName()),
                              'newick')

            distMatrix = \
                parseArray(str(obj.distanceMatrixFile.get().getFileName()))
            reordered_matrix, indices = reorderMatrix(distMatrix, tree)
            plt.figure()
            show = showMatrix(reordered_matrix, ticklabels=indices,
                              origin='upper', allticks=True)
            plt.show()

        return show

    def _viewClusterProj(self, e=None):

        obj = self.protocol
        cls = type(obj)

        if issubclass(cls, clusterPdbTrajectories):

            tree = Phylo.read(str(obj.ClusteringTree.getFileName()),
                              'newick')

            pca = loadModel(str(obj.pcaNpzFile.get().getFileName()) +
                            '.pca.npz')

            combinedTrajSet = obj.setOfTrajectories.get()
            for i, traj in enumerate(combinedTrajSet):
                if str(traj._initialPdb).endswith('.pdb'):
                    Pdb = parsePDB(str(traj._initialPdb))
                elif str(traj._initialPdb).endswith('.cif'):
                    Pdb = parseCIF(str(traj._initialPdb))
                else:
                    raise ValueError('The initial PDB should have a filename '
                                     'ending in .pdb or .cif')
                break

            protein = Pdb.select('protein and not hydrogen').copy()

            for combinedTraj in combinedTrajSet:
                combinedEns = parseDCD(str(combinedTraj.getFileName()))
                combinedEns.setCoords(protein.ca.copy())
                combinedEns.setAtoms(protein.ca)

            subgroup_color_dict = {}
            labels = []
            for group in findSubgroups(tree, float(obj.subgroupCutoff.get())):
                c = (random.random(), random.random(), random.random())
                for label in group:
                    subgroup_color_dict[label] = c
                    labels.append(label)

            labels = array(labels)

            colors = []
            for i in range(len(labels)):
                colors.append(subgroup_color_dict[labels[where(labels == str(i)
                                                               )[0][0]]])

            # plt.figure()
            show = showProjection(combinedEns, pca[:2],
                                  color=colors, markeredgewidth=0)

            ax = gca()

            projection = calcProjection(combinedEns, pca[:2])
            for n, point in enumerate(projection):
                ax.annotate(str(n), (point[0], point[1]))

            plt.show()

            return show






