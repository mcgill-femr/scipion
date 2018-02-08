
import math
import numpy as np

from pyworkflow.em import *
from prody import *
from pyworkflow.em.packages.ProDy import PRODY
from pyworkflow.utils import *
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils.path import createLink
import commands

class ProdyProt(EMProtocol):

    """ Protocol to execute functions from ProDy software"""

    _label = "ProDy protocol"

    def _defineParams(self, form):
        form.addSection(label = "Prody NMA analysis and molecular dynamics")
        form.addParam('inputStructure', PointerParam, label="Input structure",
                      important=True,
                      pointerClass='PdbFile',
                      help='The input structure can be an atomic model (true PDB) '
                           'or a pseudoatomic model(an EM volume converted into '
                           'pseudoatoms)')
        form.addParam()





