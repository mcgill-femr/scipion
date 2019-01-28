# **************************************************************************
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
# **************************************************************************

import os

import pyworkflow.tests as pwtests
from pyworkflow.tests.em.workflows import TestWorkflow
import pyworkflow.em as pwem
import pyworkflow.em.packages.relion as relion
import pyworkflow.em.packages.gctf as gctf


CPUS = os.environ.get('SCIPION_TEST_CPUS', 4)
GPUS = os.environ.get('SCIPION_TEST_GPUS', 2)


class TestWorkflowRelion3Betagal(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('relion30_tutorial')

    def _importMovies(self):
        protImport = self.newProtocol(
            pwem.ProtImportMovies,
            filesPath=self.ds.getFile('Movies/'),
            filesPattern='*.tiff',
            samplingRateMode=0,
            samplingRate=0.885,
            magnification=50000,
            scannedPixelSize=7.0,
            voltage=200,
            sphericalAberration=1.4,
            doseInitial=0.0,
            dosePerFrame=1.277,
            gainFile=self.ds.getFile("Movies/gain.mrc")
        )
        protImport.setObjLabel('import 24 movies')
        protImport.setObjComment('Relion 3 tutorial movies:\n\n'
                                 'Microscope Jeol Cryo-ARM 200\n'
                                 'Data courtesy of Takyuki Kato in the Namba '
                                 'group\n(Osaka University, Japan)')
        protImport = self.launchProtocol(protImport)

        # Validate output movies
        movies = getattr(protImport, 'outputMovies', None)
        self.assertIsNotNone(movies, "No movies were generated from the import")
        dims = movies.getDim()
        self.assertEqual((3710, 3838, 24), dims)
        self.assertEqual(24, movies.getSize())

        return protImport

    def _runRelionMc(self, protImport):
        protRelionMc = self.newProtocol(
            relion.ProtRelionMotioncor,
            objLabel='relion - motioncor',
            patchX=5, patchY=5,
            numberOfThreads=CPUS,
        )

        protRelionMc.inputMovies.set(protImport.outputMovies)
        protRelionMc = self.launchProtocol(protRelionMc)

        return protRelionMc

    def _runRelionLog(self, protRelionMc):
        protRelionLog = self.newProtocol(
            relion.ProtRelionAutopickLoG,
            objLabel='relion - autopick log',
            boxSize=250,
            minDiameter=150, maxDiameter=180,
            areParticlesWhite=False,
            numberOfThreads=CPUS,
        )

        protRelionLog.inputMicrographs.set(protRelionMc.outputMicrographs)
        protRelionLog = self.launchProtocol(protRelionLog)

        return protRelionLog

    def _runGctf(self, protMc):
        protGctf = self.newProtocol(
            gctf.ProtGctf,
            objLabel='gctf',
            lowRes=0.04, highRes=0.21,
            astigmatism=100,
            windowSize=512,
            gpuList="0"
        )

        protGctf.inputMicrographs.set(protMc.outputMicrographs)
        protGctf = self.launchProtocol(protGctf)

        return protGctf

    def _runRelionExtract(self, protPicking, protCtf):
        protRelionExtract = self.newProtocol(
            relion.ProtRelionExtractParticles,
            objLabel='relion - extract',
            boxSize=256, doRescale=True, rescaledSize=64,
            doInvert=True, doNormalize=True,
            backDiameter=50,
            numberOfMpi=CPUS/2,
            downsamplingType=0,  # Micrographs same as picking
        )

        protRelionExtract.ctfRelations.set(protCtf.outputCTF)
        protRelionExtract.inputCoordinates.set(protPicking.outputCoordinates)
        protRelionExtract = self.launchProtocol(protRelionExtract)

        return protRelionExtract

    def _runRelion2D(self, protExtract):
        protRelion2D = self.newProtocol(
            relion.ProtRelionClassify2D,
            objLabel='relion - 2d',
            maskDiameterA=200,
            numberOfClasses=100,
            extraParams='--maxsig 25',
            pooledParticles=30,
            doGpu=True,
            numberOfThreads=4,
            numberOfMpi=GPUS+1,
            allParticlesRam=True
        )

        protRelion2D.inputParticles.set(protExtract.outputParticles)
        protRelion2D = self.launchProtocol(protRelion2D)
        return protRelion2D

    def test_workflow(self):
        protImport = self._importMovies()
        protRelionMc = self._runRelionMc(protImport)
        protGctf = self._runGctf(protRelionMc)
        protRelionLog = self._runRelionLog(protRelionMc)
        protRelionExtract = self._runRelionExtract(protRelionLog, protGctf)
        protRelion2D = self._runRelion2D(protRelionExtract)


if __name__ == '__main__':
    import unittest
    unittest.main()
