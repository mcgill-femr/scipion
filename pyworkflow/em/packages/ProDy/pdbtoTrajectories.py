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
from shutil import copy

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
        form.addParam('acceptanceParam', params.FloatParam, default=0.9,
                      label='Acceptance Ratio', expertLevel =
                      pwconst.LEVEL_ADVANCED, condition = 'useFinalPdb==True',
                      help='The ratio for accepting uphill ANM-MC steps, '
                           'i.e. those that do not approach the target '
                           'structure.')
        form.addParam('maxNumSteps', params.IntParam, default=1000000,
                      label='Max number of ANM steps', expertLevel=
                      pwconst.LEVEL_ADVANCED,
                      help='Maximal number of steps in ANM-MC step.')
        form.addParam('stepCut', params.FloatParam, default=2,
                      label='RMSD cutoff for ANM-MC steps',
                      expertLevel=pwconst.LEVEL_ADVANCED,
                      help='If the total deviation from the starting '
                           'structure exceeds this value, the ANM-MC stepping '
                           'stops.')
        form.addParam('minLen', params.FloatParam, default=1,
                      label='Minimisation length  (ps)',
                      expertLevel=pwconst.LEVEL_ADVANCED,
                      help='This value is given in ps to be comparable to TMD '
                           'length. Each ps is equivalent to 500 steps.')
        form.addParam('tmdLen', params.FloatParam, default=10,
                      label='Targeted MD length (ps)',
                      expertLevel=pwconst.LEVEL_ADVANCED,
                      help='This length depends on the size of the protein '
                           'with bigger proteins needing longer. As an '
                           'indication, 10 ps was used for the leucine '
                           'transporter, which has about 500 residues.')
        '''form.addParam('spring', params.IntParam, default=20000,
                      label='Spring constant in targeted MD', expertLevel=
                      pwconst.LEVEL_ADVANCED,
                      help='In targeted molecular dynamics simulation, '
                           'the target potential is harmonic and the spring '
                           'constant term shows the force applied to a given '
                           'structure to reach the target structure')'''

        form.addSection(label='Animation')
        form.addParam('amplitude', FloatParam, default=50,
                      label="Amplitude")
        form.addParam('nframes', IntParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of frames')
        form.addParam('downsample', FloatParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      # condition=isEm
                      label='Downsample pseudoatomic structure',
                      help='Downsample factor 2 means removing one half of the '
                           'atoms or pseudoatoms.')
        form.addParam('pseudoAtomThreshold', FloatParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      # cond
                      label='Pseudoatom mass threshold',
                      help='Remove pseudoatoms whose mass is below this '
                           'threshold. '
                           'This value should be between 0 and 1.\n'
                           'A threshold of 0 implies no atom removal.')


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
                        # 'spr': self.spring.get(),
                        'cycle': self.cycleNumber.get(),
                        'minLen': self.minLen.get(),
                        'tmdLen': self.tmdLen.get(),
                        'stepCut': self.stepCut.get()
                        }
        if self.useFinalPdb.get() is True:
            self._params['finPdb'] = self.finalPdb.get().getFileName()
            self._params['cycle'] = self.defaultCycles
        else:
            self._params['finPdb'] = self._params['initPdb']

        # all_trajectories = Trajectory("all trajectories combined")
        t = time.time()
        for traj in range(self.numTrajectories.get()):
            print("Calculating trajectory %d..." %(traj+1))
            sys.stdout.flush()

            args = ('-args ' + " " + str(os.path.abspath(self._getExtraPath()))
                    + " " + "results{0}".format(str(traj+1))
                    + " " + self._params['initPdb']
                    + " " + self._params['finPdb']
                    + " " + str(self._params['cycle'])
                    + " " + str(self._params['maxDev'])
                    + " " + str(self._params['acceptParam'])
                    + " " + str(self._params['stepCut'])
                    + " " + str(self._params['minLen'])
                    + " " + str(self._params['tmdLen'])
                    + " " + str(self._params['cutoff'])
                    + " " + str(self._params['maxSteps'])
                    + " & ")

            self.runJob("VMDARGS='text with blanks' vmd -dispdev text -e " +
                        os.path.abspath(os.environ['SCIPION_HOME'] +
                        '/software/em/prody/vmd-1.9.3/lib/plugins/noarch/tcl'
                        '/comd/comd.tcl'),args)

            outputFn = self._getExtraPath('results{0}.log'.format(traj+1))
            finished = 0
            while not finished:
                time.sleep(10)
                file = open(outputFn, 'r')
                lines = file.readlines()
                for line in lines:

                    if line.find("ERROR") != -1:
                        raise OSError(line)

                    elif line.find("FINISHED") != -1:
                        finished = 1
                        print("Trajectory done.")
                        sys.stdout.flush()

                file.close()

            if self.useFinalPdb.get():
                w1_start = parsePDB(self._getExtraPath('walker1_ionized.pdb'))
                w1_traj = parseDCD(self._getExtraPath('walker1_trajectory.dcd'))
                w1_traj.setCoords(w1_start)
                w1_traj.setAtoms(w1_start.select('protein and not hydrogen'))
                writeDCD(self._getExtraPath(
                         'walker1_trajectory{:02d}_protein.dcd'.format(traj+1)),
                         w1_traj)

                w2_start = parsePDB( self._getExtraPath('walker2_ionized.pdb'))
                w2_traj = parseDCD( self._getExtraPath('walker2_trajectory.dcd'))
                w2_traj.setCoords(w2_start)
                w2_traj.setAtoms(w2_start.select('protein and not hydrogen'))
                writeDCD(self._getExtraPath(
                         'walker2_trajectory{:02d}_protein.dcd'.format(traj+1)),
                         w2_traj)

                combined_traj = parseDCD(self._getExtraPath(
                                         'walker1_trajectory{:02d}_protein.dcd'
                                         .format(traj+1)))

                for i in reversed(range(len(w2_traj))):
                    combined_traj.addCoordset(w2_traj.getConformation(i))

                writeDCD( self._getExtraPath('trajectory{:02d}.dcd'.format(traj+1)),
                          combined_traj)
                os.system('mv walker2_trajectory.dcd walker2_trajectory{:02d}'
                          '.dcd'.format(traj + 1))
            else:
                w1_start = parsePDB( self._getExtraPath('walker1_ionized.pdb'))
                w1_traj = parseDCD( self._getExtraPath('walker1_trajectory.dcd'))
                w1_traj.setCoords(w1_start)
                w1_traj.setAtoms(w1_start.select('protein and not hydrogen'))
                writeDCD(self._getExtraPath('trajectory{:02d}.dcd'.format(traj+1)),
                         w1_traj)

            # all_trajectories.addFile(self._getExtraPath('trajectory{:02d}.dcd'.
            #                                             format(traj+1)))
            os.system('mv walker1_trajectory.dcd walker1_trajectory{:02d}.dcd'.
                      format(traj+1))

            os.system('mv rmsd.txt trajectory{:02d}_rmsd.txt'.format(traj + 1))


        elapsed = time.time()-t
        print elapsed
        # writeDCD(self._getExtraPath("combined_all_trajectories.dcd"),
        #                        all_trajectories)

        '''self._insertFunctionStep('animateModesStep', n,
                                 self.amplitude.get(), self.nframes.get(),
                                 self.downsample.get(),
                                 self.pseudoAtomThreshold.get(),
                                 self.pseudoAtomRadius)'''

        self._insertFunctionStep('createOutputStep')

    '''def animateModesStep(self, numberOfModes, amplitude, nFrames, downsample,
                             pseudoAtomThreshold, pseudoAtomRadius):
            makePath(self._getExtraPath('animations'))
            self._enterWorkingDir()

            if self.structureEM:
                fn = "pseudoatoms.pdb"
                self.runJob("nma_animate_pseudoatoms.py",
                            "%s extra/vec_ani.pkl 7 %d "
                            "%f extra/animations/"
                            "animated_mode %d %d %f" % \
                            (fn, numberOfModes, amplitude, nFrames, downsample,
                             pseudoAtomThreshold), env=getNMAEnviron())
            else:
                fn = "atoms.pdb"
                self.runJob("nma_animate_atoms.py", "%s extra/vec_ani.pkl 7 %d %f "
                                                    "extra/animations/animated_mode "
                                                    "%d" % \
                            (fn, numberOfModes, amplitude, nFrames),
                            env=getNMAEnviron())

            for mode in range(7, numberOfModes + 1):
                fnAnimation = join("extra", "animations", "animated_mode_%03d"
                                   % mode)
                fhCmd = open(fnAnimation + ".vmd", 'w')
                fhCmd.write("mol new %s.pdb\n" % self._getPath(fnAnimation))
                fhCmd.write("animate style Loop\n")
                fhCmd.write("display projection Orthographic\n")
                if self.structureEM:
                    fhCmd.write("mol modcolor 0 0 Beta\n")
                    fhCmd.write("mol modstyle 0 0 Beads %f 8.000000\n"
                                % (pseudoAtomRadius))
                else:
                    fhCmd.write("mol modcolor 0 0 Index\n")
                    # fhCmd.write("mol modstyle 0 0 Beads 1.000000 8.000000\n")
                    fhCmd.write("mol modstyle 0 0 NewRibbons 1.800000 6.000000 "
                                "2.600000 0\n")
                fhCmd.write("animate speed 0.5\n")
                fhCmd.write("animate forward\n")
                fhCmd.close();

            self._leaveWorkingDir()'''

    def createOutputStep(self):
        print("Starting createOutputStep")
        sys.stdout.flush()

        setOfPDBs = self._createSetOfPDBs()
        setOfTrajectories = self._createSetOfTrajectories()
        for n in range(self.numTrajectories.get()):
            fnDcd = self._getExtraPath('trajectory{:02d}.dcd'.format(n+1))
            ens = parseDCD(fnDcd)
            atoms = parsePDB(self._getExtraPath('walker1_ionized.pdb'))
            protein = atoms.select('protein and not hydrogen').copy()
            ens.setCoords(protein)
            ens.setAtoms(protein)
            fnPdb = []
            for i, conformation in enumerate(ens):
                fnPdb.append(self._getExtraPath('trajectory{:02d}_pdb{:02d}.pdb'.format(n + 1,i + 1)))
                writePDB(fnPdb[i], ens.getConformation(i))
                pdb = PdbFile(fnPdb[i])
                setOfPDBs.append(pdb)

            traj = TrajectoryDcd(fnDcd,
                   self._getExtraPath('trajectory{:02d}_pdb01.pdb'.format(n+1)))
            setOfTrajectories.append(traj)

        self._defineOutputs(outputPDBs=setOfPDBs)
        self._defineSourceRelation(self.initialPdb.get(), setOfPDBs)

        self._defineOutputs(outputTrajs=setOfTrajectories)
        self._defineSourceRelation(self.initialPdb.get(), setOfTrajectories)

        print("Finished createOutputStep")
        sys.stdout.flush()












