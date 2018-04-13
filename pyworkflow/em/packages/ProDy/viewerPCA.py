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
            inputType = 0

            if setOfTraj is None:
                setOfTraj = obj.setOfPDBs.get()
                inputType = 1


            if inputType == 0:

                for traj in setOfTraj:
                    if str(traj._initialPdb).endswith('.pdb'):
                        Pdb = parsePDB(str(traj._initialPdb))
                    elif str(traj._initialPdb).endswith('.cif'):
                        Pdb = parseCIF(str(traj._initialPdb))
                    else:
                        raise ValueError('The initial PDB should have a filename '
                                         'ending in .pdb or .cif')
                    if traj._pseudoatoms:
                        pseudoatoms = True
                    else:
                        pseudoatoms = False
                    break

            else:
                for pdb in setOfTraj:
                    if str(pdb.getFileName()).endswith('.pdb'):
                        Pdb = parsePDB(str(pdb.getFileName()))
                    elif str(pdb.getFileName()).endswith('.cif'):
                        Pdb = parseCIF(str(pdb.getFileName()))
                    else:
                        raise ValueError(
                            'The initial PDB should have a filename '
                            'ending in .pdb or .cif')
                    break

            if pseudoatoms:
                protein = Pdb.copy()
            else:
                protein = Pdb.select('protein and not hydrogen').copy()
            combinedTrajSet = obj.combinedTrajectory
            for combinedTraj in combinedTrajSet:
                combinedEns = parseDCD(str(combinedTraj.getFileName()))
            if pseudoatoms:
                proteinca = protein
            else:
                proteinca = protein.ca
            combinedEns.setCoords(proteinca.copy())
            combinedEns.setAtoms(proteinca)

            colors = []
            for i, traj in enumerate(setOfTraj):

                c = (random.random(), random.random(), random.random())

                if inputType == 0:
                    ens = parseDCD(str(traj._filename))
                    ens.setCoords(protein)
                    ens.setAtoms(proteinca)

                    for j in range(len(ens)):
                        if j == 0:
                            colors.append((1, 0, 0))
                        else:
                            colors.append(c)
                else:
                    if i > 0:
                        colors.append(c)

            plt.figure()
            show = showProjection(combinedEns, pca[:2],
                                  color=colors, markeredgewidth=0,
                                  norm=False)

            ax = gca()

            projection = calcProjection(combinedEns, pca[:2],
                                        norm=False)

            for n, point in enumerate(projection):
                ax.annotate(str(n), (point[0], point[1]))

            '''distanceMatrix = parseArray(str(obj.distanceMatrix
                                            .getFileName()), '\t')
            plt.figure()
            showMatrix(distanceMatrix, origin='upper')'''
            plt.show()
