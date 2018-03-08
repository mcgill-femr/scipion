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

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from pyworkflow.gui.project import ProjectWindow
from pdbtoTrajectories import computePdbTrajectories
from prody import *
from pyworkflow.em.packages.xmipp3.nma import viewer_nma

OBJCMD_NMA_VMD = "Display VMD animation"

class ProdyTrajectoriesViewer(viewer_nma):

    _targets = [computePdbTrajectories]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
        self.viewer_nma._defineParams(form)

    def _getVisualizeDict(self):
        return {'displayModes': self._viewParam,
                'displayVmd': self._viewSingleMode,
                }

    def _viewSingleModes(self):
        self.viewer_nma._viewSingleModes(self._getVisualizeDict())


ProjectWindow.registerObjectCommand(OBJCMD_NMA_VMD,
                                    viewer_nma.showVmdView)