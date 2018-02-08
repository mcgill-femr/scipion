# coding: latin-1
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
Bibtex string file for Xmipp package.
"""

_bibtexStr = """


@article{Kurkcuoglu2016,
title = "ClustENM: ENM-Based Sampling of Essential Conformational Space at Full Atomic Resolution",
journal = "J. Chemical Theory and Computation",
volume = "12",
issue = "9",
pages = "4549 - 4562",
year = "2016",
doi = "10.1021/acs.jctc.6b00319",
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5088496/",
author = "Kurkcuoglu, Z. and Bahar, I. and Doruker, P.",
}
"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)