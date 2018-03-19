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

from pyworkflow.em.protocol.protocol import EMProtocol
import pyworkflow.protocol as pwprot

class joinPDBs(EMProtocol):
    _label = "Join PDBs ProDy"

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputPDBs', pwprot.params.MultiPointerParam,
                      label="Input PDBs", important=True,
                      pointerClass='SetOfPDBs,PdbFile',
                      minNumObjects=2, maxNumObjects=0,
                      help='Select two or more PDBs to be united into a set.')

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
        self.outputSet = self._createSetOfPDBs()
        for Pdb in self.inputPDBs:
            self.outputSet.append(Pdb.get())

        self._defineOutputs(SetOfPDBs=self.outputSet)
