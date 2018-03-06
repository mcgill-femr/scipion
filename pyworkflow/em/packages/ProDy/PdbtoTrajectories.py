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

SCIPION_HOME = os.environ['SCIPION_HOME']

class computePdbTrajectories(EMProtocol):
    """ Protocol to execute functions from ProDy software"""

    _label = "Generate Trajectories ProDy"

    def _defineParams(self, form):
        self.defaultCycles = 5

        form.addSection(label="Prody Trajectories")
        form.addParam('initialPdb', PointerParam, pointerClass='PdbFile',
                      label='Initial Pdb', important=False,
                      help='Choose an initial conformation using a Pdb to '
                           'compute its N trajectories')
        form.addParam('useFinalPdb', params.BooleanParam, default=False,
                      label="Use final conformation")
        form.addParam('finalPdb', PointerParam, pointerClass='PdbFile',
                      label='Final Pdb', important=True,
                      condition = 'useFinalPdb == True',
                      help='Choose a final conformation using a Pdb to '
                           'compute its N trajectories')
        form.addParam('cycleNumber', params.IntParam,
                      default=self.defaultCycles,
                      label='Number of cycles',
                      condition = 'useFinalPdb == False',
                      important=True, help='Number of cycles to cover a '
                                           'landscape')
        form.addParam('numTrajectories', params.IntParam, default=1,
                      label='Number of trajectories',
                      important=True, help='Number of trajectories to generate')
        #form.addParam('inputTopology', params.MultiPointerParam,
        #              pointerClass='SetOfVolumes,Volume',
        #              label="Input topology files", important=True,
        #              help='Select topology files')
        #form.addParam('inputParameters', params.MultiPointerParam,
        #              label="Input parameters files", important=True,
        #              help='Select parameters files')
        form.addParam('anmCutoff', params.IntParam, default=15,
                      label='ANM Cutoff', expertLevel = pwconst.LEVEL_ADVANCED,
                      help='Maximum distance that two residues are in contact.')
        form.addParam('maxDeviation', params.FloatParam, default=1.5,
                      label='Maximum deviation per normal mode (A)',
                      expertLevel = pwconst.LEVEL_ADVANCED,
                      help='The maximum deviation per step when disturbing the '
                           'protein structure in ANM-MC.')
        form.addParam('acceptanceParam', params.FloatParam, default=0.1,
                      label='Acceptance Parameter', expertLevel =
                      pwconst.LEVEL_ADVANCED, condition = 'useFinalPdb==True',
                      help='starting value for the acceptance parameter in '
                           'ANM-MC steps')
        form.addParam('maxNumSteps', params.IntParam, default=1000000,
                      label='Max number of ANM steps', expertLevel=
                      pwconst.LEVEL_ADVANCED,
                      help='Maximal number of steps in ANM-MC step.')
        form.addParam('spring', params.IntParam, default=20000,
                      label='Spring constant in targeted MD', expertLevel=
                      pwconst.LEVEL_ADVANCED,
                      help='In targeted molecular dynamics simulation, '
                           'the target potential is harmonic and the spring '
                           'constant term shows the force applied to a given '
                           'structure to reach the target structure')


    def _insertAllSteps(self):
        self._insertFunctionStep('createTrajectories')
        self._insertFunctionStep('createOutputStep')

    def createTrajectories(self):
        

        self._params = {'initPdb': self.initialPdb.get().getFileName(),
                        'numTraj': self.numTrajectories.get(),
                        'cutoff': self.anmCutoff.get(),
                        'maxDev': self.maxDeviation.get(),
                        'acceptParam': self.acceptanceParam.get(),
                        'maxSteps': self.maxNumSteps.get(),
                        'spr': self.spring.get(),
                        'cycle': self.cycleNumber.get()
                        }

        if self.useFinalPdb.get() is True:
            self._params['finPdb'] = self.finalPdb.get().getFileName()
            self._params['cycle'] = self.defaultCycles
        else:
            self._params['finPdb'] = self._params['initPdb']

        for traj in range(self.numTrajectories.get()):

            os.system("env VMDARGS='text with blanks' vmd -dispdev text -e " +
                            SCIPION_HOME + "/software/em/prody/vmd-1.9.3/lib/"
                            "plugins/noarch/tcl/comd/comd.tcl -args "
                            + os.path.abspath(self._getExtraPath()) + " "
                            + "results{0}".format(str(traj+1))
                            + " " + self._params['initPdb']
                            + " " + self._params['finPdb']
                            + " " + str(self._params['cycle'])
                            + " " + str(self._params['maxDev'])
                            + " 2>> " + self._getLogsPath('run.stderr')
                            #+ " 1>> " + self._getLogsPath('run.stdout')
                            + " & ")

            print("Calculating trajectory %d..." %(traj+1))
            sys.stdout.flush()
            outputFn = self._getExtraPath('results{0}.log'.format(traj+1))
            finished = 0
            while not finished:
                time.sleep(20)
                file = open(outputFn, 'r')
                lines = file.readlines()
                for line in lines:

                    if line.find("ERROR") != -1:
                        raise OSError('The job had an error and'
                                      'could not calculate the '
                                      'trajectory. Please see '
                                      'results{0}.log'.format(traj+1))

                    elif line.find("FINISHED") != -1:
                        finished = 1
                        print("Trajectory done.")
                        sys.stdout.flush()

                file.close()

            if self.useFinalPdb.get():
                w1_start = parsePDB( self._getExtraPath('walker1_ionized.pdb'))
                w1_traj = parseDCD( self._getExtraPath('walker1_trajectory.dcd'))
                w1_traj.setCoords(w1_start)
                w1_traj.setAtoms(w1_start.select('protein and not hydrogen'))
                writeDCD( self._getExtraPath(
                    'walker1_trajectory_protein.dcd'), w1_traj)

                w2_start = parsePDB( self._getExtraPath('walker2_ionized.pdb'))
                w2_traj = parseDCD( self._getExtraPath('walker2_trajectory.dcd'))
                w2_traj.setCoords(w2_start)
                w2_traj.setAtoms(w2_start.select('protein and not hydrogen'))
                writeDCD( self._getExtraPath(
                    'walker2_trajectory_protein.dcd'), w2_traj)

                combined_traj = parseDCD( self._getExtraPath('walker1_trajectory_protein.dcd'))

                for i in reversed(range(len(w2_traj))):
                    combined_traj.addCoordset(w2_traj.getConformation(i))

                writeDCD( self._getExtraPath('trajectory.dcd'), combined_traj)
            else:
                w1_start = parsePDB( self._getExtraPath('walker1_ionized.pdb'))
                w1_traj = parseDCD( self._getExtraPath('walker1_trajectory.dcd'))
                w1_traj.setCoords(w1_start)
                w1_traj.setAtoms(w1_start.select('protein and not hydrogen'))
                writeDCD( self._getExtraPath('trajectory.dcd'), w1_traj)


            os.system("mv  " + self._getExtraPath("trajectory.dcd") + " " +
                self._getExtraPath("trajectory%i.dcd" %(traj+1)))

    def createOutputStep(self):
        print("En createOutputStep")
        sys.stdout.flush()












