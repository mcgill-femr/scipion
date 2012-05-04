#!/usr/bin/env xmipp_python
'''
#/***************************************************************************
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************
'''

from xmipp import MDL_CTF_SAMPLING_RATE, MDL_CTF_VOLTAGE, MDL_CTF_DEFOCUSU, MDL_CTF_DEFOCUSV, \
MDL_CTF_DEFOCUS_ANGLE, MDL_CTF_CS, MDL_CTF_CA, MDL_CTF_Q0, MDL_CTF_K, label2Str, MetaData


def convertCtfparam(oldCtf):
    '''Convert the old format (Xmipp2.4) of the CTF 
    and return a new MetaData'''
    conversionDict = {'sampling_rate': MDL_CTF_SAMPLING_RATE, 
                      'voltage':       MDL_CTF_VOLTAGE, 
                      'defocusU':      MDL_CTF_DEFOCUSU, 
                      'defocusV':      MDL_CTF_DEFOCUSV, 
                      'azimuthal_angle':      MDL_CTF_DEFOCUS_ANGLE, 
                      'spherical_aberration': MDL_CTF_CS,
                      'chromatic_aberration': MDL_CTF_CA,
                      'Q0':                   MDL_CTF_Q0, 
                      'K':                    MDL_CTF_K}
    
    f = open(oldCtf)
    md = MetaData()
    md.setColumnFormat(False)
    objId = md.addObject()
    
    for line in f:
        line = line.strip()
        if len(line) and not line.startswith('#'):
            parts = line.split('=')
            old_key = parts[0].strip()
            value = float(parts[1].strip())
            label = conversionDict[old_key]
            md.setValue(label, value, objId)

    f.close()
  
    #Change sign of defocusU, defocusV and Q0
    labels = [label2Str(l) for l in [MDL_CTF_DEFOCUSU, MDL_CTF_DEFOCUSV, MDL_CTF_Q0]]
    exps = ["%s=%s*-1" % (l, l) for l in labels]
    expression = ','.join(exps)
    
    md.operate(expression)
    return md
    
    
    