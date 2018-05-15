# **************************************************************************
# *
# * Authors:     Javier Vargas Balbuena (javier.vargasbalbuena@mcgill.ca)
# *
# * Department of Anatomy and Cell Biology, McGill University
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


from pyworkflow.protocol.params import (PointerParam, FloatParam,  
                                        StringParam, BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
from pyworkflow.em.data import Volume


class Prot3DVariance(ProtAnalysis3D):
    """    
    Obtains the 3D variance in real and fourier spaces from a given set of particles.
    First obtained a 3D reconstruction using Xmipp reconstruct Fourier method.
    """
    _label = '3d_variance'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles', pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help='Select the input images from the project.')     
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        form.addParam('maxRes', FloatParam, default=-1,
                      label="Maximum resolution (A)",  
                      help='Maximum resolution (in Angstrom) to consider \n'
                           'in Fourier space (default Nyquist).\n'
                           'Param *--maxres* in Xmipp.')
        form.addParam('doCTFCorrection', BooleanParam, default=False,  important=False,
                      label='Use CTF information?',
                      help="The images will be deconvolved by the CTF information.")
        form.addParam('Mask', PointerParam, pointerClass='VolumeMask',allowsNull=True, 
                      label="Binary Mask",
                      help='Binary mask that determines which points are specimen'
                      ' and which are not')  
        form.addParam('pad', FloatParam, default=2,
                      label="Padding factor")
        form.addParam('extraParams', StringParam, default='', expertLevel=LEVEL_ADVANCED,
                      label='Extra parameters: ', 
                      help='Extra parameters to *xmipp_reconstruct_fourier-xmipp_volume_variability* programs:\n')

        form.addParallelSection(threads=1, mpi=4)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_xmd': self._getExtraPath('input_particles.xmd'),
            'output_volume': self._getPath('output_volume.vol')
            }
        self._updateFilenamesDict(myDict)

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        self._insertReconstructStep()
        self._insertVarianceStep()
        self._insertFunctionStep('createOutputStep')
        
    def _insertReconstructStep(self):
        #imgSet = self.inputParticles.get()

        self.params =  '  -i %s' % self._getFileName('input_xmd')
        self.params += '  -o %s' % self._getFileName('output_volume')
        self.params += ' --sym %s' % self.symmetryGroup.get()
        maxRes = self.maxRes.get()
        if maxRes == -1:
            digRes = 0.5
        else:
            digRes = self.inputParticles.get().getSamplingRate() / self.maxRes.get()
        
        if self.doCTFCorrection:
            self.params +=  '  --useCTF'                  
            
        self.params += ' --max_resolution %0.3f' %digRes
        self.params += ' --padding %0.3f' % self.pad.get()
        self.params += ' --thr %d' % self.numberOfThreads.get()
        self.params += ' --sampling %f' % self.inputParticles.get().getSamplingRate()
        self.params += ' %s' % self.extraParams.get()
        
        self._insertFunctionStep('reconstructStep', self.params)

    def _insertVarianceStep(self):
        
        if (self.Mask.get() is not None):
            self.params += '   --mask %s' %  self.Mask.get().getFileName()
                   
        self._insertFunctionStep('varianceStep', self.params)
        
    #--------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        particlesMd = self._getFileName('input_xmd')
        imgSet = self.inputParticles.get()
        #TODO: This only writes metadata what about binary file
        #it should
        writeSetOfParticles(imgSet, particlesMd)

    def reconstructStep(self, params):
        """ Create the input file in STAR format as expected by Xmipp.
        If the input particles comes from Xmipp, just link the file. 
        """
        self.runJob('xmipp_reconstruct_fourier', params)

    def varianceStep(self, params):
        self.runJob('xmipp_volume_variability', params)
            
    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        volume_reconst = Volume()
        volume_reconst.setFileName(self._getFileName('output_volume'))
        volume_reconst.setSamplingRate(imgSet.getSamplingRate())

        volume_variance = Volume()
        volume_variance.setFileName(self._getFileName('variability_variance'))
        volume_variance.setSamplingRate(imgSet.getSamplingRate())
        
        self._defineOutputs(outputVolume=volume_reconst)
        self._defineSourceRelation(self.inputParticles, volume_reconst)

        self._defineOutputs(outputVolume=volume_variance)
        self._defineSourceRelation(self.inputParticles, volume_variance)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    #--------------------------- UTILS functions --------------------------------------------
