# **************************************************************************
# *
# * Authors:     Amaya Jimenez (ajimenez@cnb.csic.es)
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


from pyworkflow.protocol.params import *
from pyworkflow.em.protocol import ProtImportFiles
from pyworkflow.em.data import TrajectoryDcd
from os import listdir
from os.path import isfile, join
from shutil import copy
from prody import *
import glob


class ProtImportTrajectores(ProtImportFiles):
    """ Protocol to import trajectories for Prody protocols """
    _label = 'import trajectories'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Import')
        form.addParam('filesPath', PathParam,
                      label="Files directory",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/Particles/data/day??_micrographs/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/Particles/data/day*_micrographs/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/Particles/data/day#_micrographs/\n"
                           "'#' represents one digit that will be used as "
                           "micrograph ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")
        form.addParam('filesPattern', StringParam,
                      label='Pattern',
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.")
        form.addParam('initialPdb', PointerParam, pointerClass='PdbFile',
                      label='Initial Pdb', important=True,
                      help='The initial conformation for the trajectories to '
                           'import')


    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importTrajectoriesStep', self.getPattern())
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def importTrajectoriesStep(self, pattern):
        """ Copy trajectories matching the filename pattern"""

        setTrajOut = self._createSetOfTrajectories()

        print("path", self.filesPath.get())
        print("pattern", self.filesPattern.get())

        directory = self.filesPath.get()
        pattern = self.filesPattern.get()
        fnList = glob.glob1(directory, pattern)
        print(fnList)
        fnPdb = self.initialPdb.get().getFileName()
        copy(fnPdb, self._getExtraPath('initialPdb.pdb'))

        for fn in fnList:
            wholeFn = join(directory, fn)
            if isfile(wholeFn):
                copy(wholeFn, self._getExtraPath(fn))
                traj = TrajectoryDcd(self._getExtraPath(fn),
                                     self._getExtraPath('initialPdb.pdb'))
                setTrajOut.append(traj)

        self._defineOutputs(outputTrajs=setTrajOut)

    def _summary(self):
        summary = []
        return summary    
    
    def _methods(self):
        return []
    
    def _stepsCheck(self):
        # Just to avoid the stream checking inherited from ProtCTFMicrographs
        pass

    def _validate(self):
        errors = []
        if not self.getPattern():
            errors.append("The path and pattern can not be both empty!!!")
        else:
            # Just check the number of files matching the pattern
            self.getMatchFiles()
            if self.numberOfFiles == 0:
                errors.append("There are no files matching the pattern %s"
                              % self.getPattern())

        return errors
