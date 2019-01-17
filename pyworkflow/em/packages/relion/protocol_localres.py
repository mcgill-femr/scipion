# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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

from os.path import exists
import pyworkflow as pw
from pyworkflow.protocol.params import (PointerParam, FloatParam, FileParam,
                                        IntParam, LabelParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ImageHandler

from protocol_postprocess import ProtRelionPostprocess
import convert


class ProtRelionLocalRes(ProtRelionPostprocess):
    """
    Relion local-resolution estimation protocol.

    This program basically performs a series of post-processing operations
    with a small soft, spherical mask that is moved over the entire map,
    while using phase-randomisation to estimate the convolution effects
    of that mask.
    """
    _label = 'local resolution'
    _lastUpdateVersion = pw.VERSION_1_2

    @classmethod
    def isDisabled(cls):
        return convert.isVersion1()
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
                 'half1': self._getInputPath("relion_half1_class001_unfil.mrc"),
                 'half2': self._getInputPath("relion_half2_class001_unfil.mrc"),
                 'outputVolume': self._getExtraPath('relion_locres_filtered.mrc'),
                 'resolMap': self._getExtraPath('relion_locres.mrc')
                 }

        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        
        form.addSection(label='Input')
        form.addParam('protRefine', PointerParam,
                      pointerClass="ProtRefine3D",
                      label='Select a previous refinement protocol',
                      help='Select any previous refinement protocol to get the '
                           '3D half maps. Note that it is recommended that the '
                           'refinement protocol uses a gold-standard method.')
        form.addParam('mtf', FileParam,
                       label='MTF-curve file',
                       help='User-provided STAR-file with the MTF-curve '
                            'of the detector.'
                            'Relion param: <--mtf>')
        form.addParam('bfactor', FloatParam, default=-250,
                       label='Provide B-factor:',
                       help='Probably, the overall B-factor as was '
                            'estimated in the postprocess is a useful '
                            'value for here. Use negative values for '
                            'sharpening. Be careful: if you over-sharpen '
                            'your map, you may end up interpreting '
                            'noise for signal!')

        form.addSection(label='LocalRes')
        form.addParam('Msg', LabelParam,
                      label='Select Advanced level if you want to adjust the '
                            'parameters')
        form.addParam('locResSamp', IntParam, default=25,
                      label='Sampling rate (A)',
                      expertLevel=LEVEL_ADVANCED,
                      help='Sampling rate (in Angstroms) with which to '
                           'sample the local-resolution map')
        form.addParam('locResMaskRad', IntParam, default=-1,
                      label='Mask radius (A)',
                      expertLevel=LEVEL_ADVANCED,
                      help='Radius (in A) of spherical mask for '
                           'local-resolution map (default = 0.5*sampling)')
        form.addParam('locResEdgeWidth', IntParam, default=-1,
                      label='Edge width (A)',
                      expertLevel=LEVEL_ADVANCED,
                      help='Width of soft edge (in A) on masks for '
                           'local-resolution map (default = sampling)')
        form.addParam('locResRand', FloatParam, default=25.0,
                      label='Randomize phases from (A)',
                      expertLevel=LEVEL_ADVANCED,
                      help='Randomize phases from this resolution (in A)')
        form.addParam('locResMin', IntParam, default=50,
                      label='Lowest res limit (A)',
                      expertLevel=LEVEL_ADVANCED,
                      help='Lowest local resolution allowed (in A)')

        form.addParallelSection(threads=0, mpi=1)
    
    # -------------------------- INSERT steps functions -----------------------

    # ------------------------- STEPS functions -------------------------------
    def convertInputStep(self, protId):
        pw.utils.makePath(self._getInputPath())

        protRef = self.protRefine.get()
        _, half1, half2 = protRef.getFinalVolumes()
        ih = ImageHandler()
        ih.convert(half1, self._getFileName("half1"))
        ih.convert(half2, self._getFileName("half2"))
    
    # -------------------------- INFO functions -------------------------------
    def _summary(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volume not ready yet.")
        else:
            output = self.outputVolume
            summary.append("%s: Output volume was locally filtered "
                           "and sharpened" % self.getObjectTag(output))
        return summary
    
    # -------------------------- UTILS functions ------------------------------
    def _defineParamDict(self):
        """ Define all parameters to run relion_postprocess"""
        volume = self.protRefine.get().outputVolume
        # It seems that in Relion3 now the input should be the map
        # filename and not the prefix as before
        if convert.isVersion3():
            inputFn = self._getFileName('half1')
        else:
            inputFn = self._getInputPath("relion")

        self.paramDict = {'--i': inputFn,
                          '--o': self._getExtraPath('relion'),
                          '--angpix': volume.getSamplingRate(),
                          '--adhoc_bfac': self.bfactor.get(),
                          '--locres': '',
                          # Expert options
                          '--locres_sampling': self.locResSamp.get(),
                          '--locres_maskrad': self.locResMaskRad.get(),
                          '--locres_edgwidth': self.locResEdgeWidth.get(),
                          '--locres_randomize_at': self.locResRand.get(),
                          '--locres_minres': self.locResMin.get()
                          }

        mtfFile = self.mtf.get()
        if mtfFile:
            self.paramDict['--mtf'] = mtfFile

    def _getRelionMapFn(self, fn):
        return fn.split(':')[0]
