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
from computePCATrajectory import computeModesPcaPdb
from prody import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca


class ProdyViewerPca(Viewer):
    """ Wrapper to visualize Pdb to SAXS. """
    _targets = [computeModesPcaPdb]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, **args):
        Viewer.__init__(self, **args)

    def visualize(self, obj, **args):

        cls = type(obj)

        if issubclass(cls, computeModesPcaPdb):
            pca = loadModel(str(obj.pcaNpzFile.getFileName()) + '.pca.npz')

            setOfTraj = obj.setOfTrajectories.get()

            for i, traj in enumerate(setOfTraj):
                Pdb = parsePDB(str(traj._initialPdb))
                break

            protein = Pdb.select('protein and not hydrogen').copy()

            combinedTrajSet = obj.combinedTrajectory
            for combinedTraj in combinedTrajSet:
                combinedEns = parseDCD(str(combinedTraj.getFileName()))

            combinedEns.setCoords(protein.ca.copy())
            combinedEns.setAtoms(protein.ca)

            colors = []
            for i, traj in enumerate(obj.setOfTrajectories.get()):

                c = (random.random(), random.random(), random.random())

                ens = parseDCD(str(traj._filename))
                ens.setCoords(protein)
                ens.setAtoms(protein.ca)

                for j in range(len(ens)):
                    if j == 0:
                        colors.append((1,0,0))
                    else:
                        colors.append(c)

            show = showProjection(combinedEns, pca[:2],
                                  color=colors, markeredgewidth=0)

            ax = gca()

            projection = calcProjection(combinedEns, pca[:2])
            for n, point in enumerate(projection):
                ax.annotate(str(n), (point[0], point[1]))

            distanceMatrix = parseArray(str(obj.distanceMatrix
                                            .getFileName()),'\t')
            plt.figure()
            showMatrix(distanceMatrix)
            plt.show()





