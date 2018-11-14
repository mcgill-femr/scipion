# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
# ******************************************************************************

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.em as em
import pyworkflow.em.metadata as md

import convert


class ProtRelionBayesianPolishing(em.ProtParticles):
    """
    Wrapper protocol for the Relion's Bayesian Polishing.

    As of release 3.0, Relion also implements a new Bayesian approach to beam
    induced motion correction. This approach aims to optimise a regularised
    likelihood, which allows us to associate with each hypothetical set of
    particle trajectories a prior likelihood that favors spatially coherent
    and temporally smooth motion without imposing any hard constraints.
    The smoothness prior term requires three parameters that describe the
    statistics of the observed motion. To estimate the prior that yields the
    best motion tracks for this particular dataset, we can first run the
    program in 'training mode'. Once the estimates have been obtained, one
    can then run the program again to fit tracks for the motion of all
    particles in the data set and to produce adequately weighted averages of
    the aligned movie frames.

    """

    _label = 'bayesian polishing'

    OP_TRAIN = 0
    OP_POLISH = 1

    @classmethod
    def isDisabled(cls):
        return not convert.isVersion3()

    def _defineParams(self, form):
        form.addSection(label='Input')

        # TODO: conditions on movies?
        form.addParam('inputMovies', params.PointerParam, pointerClass='SetOfMovies',
                      important=True,
                      label='Input movies',
                      help='')
        #TODO: conditions on particles?
        form.addParam('inputParticles', params.PointerParam,
                      important=True,
                      label='Input particles',
                      pointerClass='SetOfParticles',
                      help='Provide a set of particles for local CTF refinement.')
        form.addParam('inputPostprocess', params.PointerParam,
                      important=True,
                      label='Input Postprocess',
                      pointerClass='ProtRelionPostprocess',
                      help='Select a PostProcess job. The mask used for this '
                           'postprocessing will be applied to the unfiltered '
                           'half-maps and should encompass the entire complex. '
                           'The resulting FSC curve will be used for weighting '
                           'the different frequencies. ')
        line = form.addLine('Movie frames',
                             help='First and last frames to take into account '
                                  'in motion fit and combination step. '
                                  '(starts counting at 1 and 0 as last '
                                  'means util the last frame in the movie).')
        line.addParam('frame0', params.IntParam, default=1,
                      label='first')
        line.addParam('frameN', params.IntParam, default=0,
                      label='last')

        form.addSection(label='Train or Polish')
        form.addParam('operation', params.EnumParam, default=0,
                      choices=[' Train optimal parameters ',
                               ' Perform particle polishing '],
                      display=params.EnumParam.DISPLAY_COMBO,
                      label='Operation',
                      help="If *train optimal parameters* , then "
                           "relion_motion_refine will estimate optimal "
                           "parameter values for the three sigma values above "
                           "on a subset of the data (determined by the minimum "
                           "number of particles to be used below).\n\n"
                           "If *perform particle polishing* then "
                           "relion_motion_refine will be run to estimate "
                           "per-particle motion-tracks using the parameters "
                           "below, and polished particles will be generated. ")

        condTrain = "operation==%s" % self.OP_TRAIN
        group = form.addGroup('Train', condition=condTrain)
        group.addParam('fractionFourierPx', params.FloatParam, default=0.5,
                      label='Fraction of Fourier pixels for testing',
                      help="This fraction of Fourier pixels (at higher "
                           "resolution) will be used for evaluation of the "
                           "parameters (test set), whereas the rest (at lower "
                           "resolution) will be used for parameter estimation "
                           "itself (work set).")
        group.addParam('numberOfParticles', params.IntParam, default=10000,
                      label='Use this many particles',
                      help='Use at least this many particles for the '
                           'meta-parameter optimisation. The more particles '
                           'the more expensive in time and computer memory '
                           'the calculation becomes, but the better the results '
                           'may get.')

        condPolish = "operation==%s" % self.OP_POLISH
        group = form.addGroup('Polish', condition=condPolish)
        group.addParam('sigmaVel', params.FloatParam, default=0.2,
                       label='Sigma for velocity (A/dose)',
                       help='Standard deviation for the velocity regularisation. '
                            'Smaller values requires the tracks to be shorter.')
        group.addParam('sigmaDiv', params.FloatParam, default=5000,
                       label='Sigma for divergence (A)',
                       help='Standard deviation for the divergence of tracks '
                            'across the micrograph. Smaller values requires '
                            'the tracks to be spatially more uniform in a '
                            'micrograph.')
        group.addParam('sigmaAcc', params.FloatParam, default=2,
                       label='Sigma for acceleration (A/dose)',
                       help='Standard deviation for the acceleration '
                            'regularisation. Smaller values requires the '
                            'tracks to be straighter.')
        line = group.addLine("Resolution for B-factor fit (A)",
                             help='The minimum and maximum spatial frequencies '
                                  '(in Angstrom) used in the B-factor fit.'
                                  'If a negative value is given as the maximum,'
                                  'it is determined from the input FSC curve.')
        line.addParam('minResBfactor', params.FloatParam, default=20, label='min')
        line.addParam('maxResBfactor', params.FloatParam, default=-1, label='max')

        form.addParallelSection(threads=1, mpi=1)

    # -------------------------- STEPS functions -------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('trainOrPolishStep')
        #self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        inputParts = self.inputParticles.get()
        imgStar = self._getPath('input_particles.star')
        inputFolder = self._getPath('input')
        pwutils.makePath(inputFolder)

        self.info("Converting set from '%s' into '%s'" %
                  (inputParts.getFileName(), imgStar))

        convert.writeSetOfParticles(inputParts, imgStar, inputFolder,
                                    alignType=em.ALIGN_PROJ,
                                    fillMagnification=True,
                                    fillRandomSubset=True)

    def trainOrPolishStep(self):
        """
        Train:
        `which relion_motion_refine`
        --i CtfRefine/job024/particles_ctf_refine.star
        --f PostProcess/job023/postprocess.star
        --corr_mic MotionCorr/job002/corrected_micrographs.star
        --m1 Refine3D/job021/run_half1_class001_unfil.mrc
        --m2 Refine3D/job021/run_half2_class001_unfil.mrc
        --mask MaskCreate/job022/mask.mrc
        --first_frame 1 --last_frame -1
        --o Polish/job030/

        --min_p 5000 --eval_frac 0.5 --align_frac 0.5 --params3  --j 16

        Polish:
        `which relion_motion_refine`
        --i CtfRefine/job024/particles_ctf_refine.star
        --f PostProcess/job023/postprocess.star
        --corr_mic MotionCorr/job002/corrected_micrographs.star
        --m1 Refine3D/job021/run_half1_class001_unfil.mrc
        --m2 Refine3D/job021/run_half2_class001_unfil.mrc
        --mask MaskCreate/job022/mask.mrc
        --first_frame 1 --last_frame -1
        --o Polish/job030/

        --s_vel 0.2 --s_div 5000 --s_acc 2 --combine_frames --bfac_minfreq 20 --bfac_maxfreq -1 --j 16 bbbb

        """
        args = "--i %s " % self._getPath('input_particles.star')
        args += "--o %s " % self._getExtraPath()
        postStar = self.inputPostprocess.get()._getExtraPath('postprocess.star')
        args += "--f %s " % postStar
        args += "--m1 %s --m2 %s --mask %s " % convert.getVolumesFromPostprocess(postStar)
        args += "--corr_mic %s " % self._getInputPath("corrected_micrographs.star")
        args += "--first_frame %d --last_frame %d " % (self.frame0, self.frameN)

        if self.operation == self.OP_TRAIN:
            args += "--min_p %d " % self.numberOfParticles
            args += "--eval_frac %0.3f " % self.fractionFourierPx
            args += "--align_frac %0.3f " % self.fractionFourierPx
            args += "--params3 "
        else:  # OP_POLISH
            args += "--s_vel %0.3f " % self.sigmaVel
            args += "--s_div %0.3f " % self.sigmaDiv
            args += "--s_acc %0.3f " % self.sigmaAcc
            args += "--combine_frames "
            args += ""

        if self.doCtfFitting:
            args += "--fit_defocus "

        fitAstig = self.fitAstig.get()
        if fitAstig == 1:
            args += "--glob_astig "
        elif fitAstig == 2:
            args += "--astig "

        if self.fitMicPhaseShift:
            args += "--fit_phase "

        if self.doBeamtiltEstimation:
            args += "--fit_beamtilt "

        args += "--j %d " % self.numberOfThreads

        self.info("runJob: program: %s, args: %s" % ("relion_motion_refine", args))
        #self.runJob("relion_motion_refine", args)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)

        outImgsFn = self._getExtraPath('particles_ctf_refine.star')
        imgSet.setAlignmentProj()
        rowIterator = md.iterRows(outImgsFn, sortByLabel=md.RLN_IMAGE_ID)
        outImgSet.copyItems(imgSet,
                            updateItemCallback=self._updateItemCtf,
                            itemDataIterator=rowIterator)

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)

    def _updateItemCtf(self, particle, row):
        particle.setCTF(convert.rowToCtfModel(row))
        #TODO: Add other field from the .star file when other options?

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        errors = []
        return errors

    def _getInputPath(self, *paths):
        return self._getPath('input', *paths)

