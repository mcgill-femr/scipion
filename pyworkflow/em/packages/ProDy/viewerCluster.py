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


class ProdyViewerCluster(Viewer):
    """ Wrapper to visualize Pdb to SAXS. """
    _targets = [clusterPdbTrajectories]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, **args):
        Viewer.__init__(self, **args)

    def visualize(self, obj, **args):

        cls = type(obj)

        if issubclass(cls, clusterPdbTrajectories):

            tree = Phylo.read(str(obj.ClusteringTree.getFileName()),
                              'newick')

            distMatrix = parseArray(str(obj.distanceMatrixFile.get().getFileName()))

            subgroup_color_dict = {}
            labels = []
            for group in findSubgroups(tree, float(obj.subgroupCutoff.get())):
                c = (random.random(), random.random(), random.random())
                for label in group:
                    subgroup_color_dict[label] = c
                    labels.append(label)

            labels = array(labels)

            showTree(tree, format='plt', label_colors=subgroup_color_dict)

            plt.figure()
            reordered_matrix, indices = reorderMatrix(distMatrix, tree)

            if len(indices) < 15:
                all_ticks = True
            else:
                all_ticks = False

            showMatrix(reordered_matrix, ticklabels=indices, origin='upper',
                       allticks=all_ticks)

            pca = loadModel(str(obj.pcaNpzFile.get().getFileName()) +
                            '.pca.npz')

            combinedTrajSet = obj.setOfTrajectories.get()
            for i, traj in enumerate(combinedTrajSet):
                Pdb = parsePDB(str(traj._initialPdb))
                break

            protein = Pdb.select('protein and not hydrogen').copy()

            for combinedTraj in combinedTrajSet:
                combinedEns = parseDCD(str(combinedTraj.getFileName()))
                combinedEns.setCoords(protein.ca.copy())
                combinedEns.setAtoms(protein.ca)

            colors = []
            for i in range(len(labels)):
                colors.append(subgroup_color_dict[labels[where(labels == str(i)
                                                               )[0][0]]])

            plt.figure()
            show = showProjection(combinedEns, pca[:2],
                                  color=colors, markeredgewidth=0,
                                  norm=False)

            ax = gca()

            projection = calcProjection(combinedEns, pca[:2],
                                        norm=False)
            for n, point in enumerate(projection):
                ax.annotate(str(n), (point[0], point[1]))

            plt.show()





