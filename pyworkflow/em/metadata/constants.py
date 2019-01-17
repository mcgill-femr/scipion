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
Add constants from xmipp module.
"""
import xmipp

#------------------ Xmipp METADATA LABELS  --------------------------------------

AGGR_COUNT = xmipp.AGGR_COUNT
AGGR_MAX = xmipp.AGGR_MAX
AGGR_SUM = xmipp.AGGR_SUM
AGGR_AVG = xmipp.AGGR_AVG
AGGR_MIN = xmipp.AGGR_MIN

UNION = xmipp.UNION
UNION_DISTINCT = xmipp.UNION_DISTINCT
INTERSECTION = xmipp.INTERSECTION
SUBSTRACTION = xmipp.SUBSTRACTION
INNER_JOIN = xmipp.INNER_JOIN
LEFT_JOIN = xmipp.LEFT_JOIN
NATURAL_JOIN = xmipp.NATURAL_JOIN
OUTER_JOIN = xmipp.OUTER_JOIN
INNER = xmipp.INNER
LEFT = xmipp.LEFT
OUTER = xmipp.OUTER
NATURAL = xmipp.NATURAL
EQ = xmipp.EQ
NE = xmipp.NE
GT = xmipp.GT
LT = xmipp.LT
GE = xmipp.GE
LE = xmipp.LE
MD_OVERWRITE = xmipp.MD_OVERWRITE
MD_APPEND = xmipp.MD_APPEND

# ----- Xmipp Metadata Labels
MDL_UNDEFINED = xmipp.MDL_UNDEFINED
MDL_FIRST_LABEL = xmipp.MDL_FIRST_LABEL
MDL_OBJID = xmipp.MDL_OBJID
MDL_ANGLE_PSI2 = xmipp.MDL_ANGLE_PSI2
MDL_ANGLE_PSI = xmipp.MDL_ANGLE_PSI
MDL_ANGLE_PSI_DIFF = xmipp.MDL_ANGLE_PSI_DIFF
MDL_ANGLE_ROT2 = xmipp.MDL_ANGLE_ROT2
MDL_ANGLE_ROT = xmipp.MDL_ANGLE_ROT
MDL_ANGLE_ROT_DIFF = xmipp.MDL_ANGLE_ROT_DIFF
MDL_ANGLE_TILT2 = xmipp.MDL_ANGLE_TILT2
MDL_ANGLE_TILT = xmipp.MDL_ANGLE_TILT
MDL_ANGLE_TILT_DIFF = xmipp.MDL_ANGLE_TILT_DIFF
MDL_ANGLE_DIFF = xmipp.MDL_ANGLE_DIFF
MDL_ANGLE_Y = xmipp.MDL_ANGLE_Y
MDL_ANGLE_Y2 = xmipp.MDL_ANGLE_Y2
MDL_AVG = xmipp.MDL_AVG
MDL_AVG_CHANGES_ORIENTATIONS = xmipp.MDL_AVG_CHANGES_ORIENTATIONS
MDL_AVG_CHANGES_OFFSETS = xmipp.MDL_AVG_CHANGES_OFFSETS
MDL_AVG_CHANGES_CLASSES = xmipp.MDL_AVG_CHANGES_CLASSES

MDL_BGMEAN = xmipp.MDL_BGMEAN
MDL_BLOCK_NUMBER = xmipp.MDL_BLOCK_NUMBER

MDL_CLASS_COUNT = xmipp.MDL_CLASS_COUNT
MDL_CLASS_PERCENTAGE = xmipp.MDL_CLASS_PERCENTAGE
MDL_CLASSIFICATION_DATA = xmipp.MDL_CLASSIFICATION_DATA
MDL_CLASSIFICATION_DATA_SIZE = xmipp.MDL_CLASSIFICATION_DATA_SIZE
MDL_CLASSIFICATION_DPR_05 = xmipp.MDL_CLASSIFICATION_DPR_05
MDL_CLASSIFICATION_FRC_05 = xmipp.MDL_CLASSIFICATION_FRC_05
MDL_CLASSIFICATION_INTRACLASS_DISTANCE = xmipp.MDL_CLASSIFICATION_INTRACLASS_DISTANCE

MDL_COLOR = xmipp.MDL_COLOR
MDL_COMMENT = xmipp.MDL_COMMENT
MDL_CONTINUOUS_FLIP = xmipp.MDL_CONTINUOUS_FLIP
MDL_CONTINUOUS_GRAY_A = xmipp.MDL_CONTINUOUS_GRAY_A
MDL_CONTINUOUS_GRAY_B = xmipp.MDL_CONTINUOUS_GRAY_B
MDL_CONTINUOUS_X = xmipp.MDL_CONTINUOUS_X
MDL_CONTINUOUS_Y = xmipp.MDL_CONTINUOUS_Y
MDL_COST = xmipp.MDL_COST
MDL_COUNT = xmipp.MDL_COUNT
MDL_COUNT2 = xmipp.MDL_COUNT2
MDL_CRYSTAL_LATTICE_A = xmipp.MDL_CRYSTAL_LATTICE_A
MDL_CRYSTAL_LATTICE_B = xmipp.MDL_CRYSTAL_LATTICE_B
MDL_CRYSTAL_DISAPPEAR_THRE = xmipp.MDL_CRYSTAL_DISAPPEAR_THRE
MDL_CRYSTAL_SHFILE = xmipp.MDL_CRYSTAL_SHFILE
MDL_CRYSTAL_ORTHO_PRJ = xmipp.MDL_CRYSTAL_ORTHO_PRJ
MDL_CRYSTAL_PROJ = xmipp.MDL_CRYSTAL_PROJ
MDL_CRYSTAL_CELLX = xmipp.MDL_CRYSTAL_CELLX
MDL_CRYSTAL_CELLY = xmipp.MDL_CRYSTAL_CELLY
MDL_CRYSTAL_SHIFTX = xmipp.MDL_CRYSTAL_SHIFTX
MDL_CRYSTAL_SHIFTY = xmipp.MDL_CRYSTAL_SHIFTY
MDL_CRYSTAL_SHIFTZ = xmipp.MDL_CRYSTAL_SHIFTZ
MDL_CTF_INPUTPARAMS = xmipp.MDL_CTF_INPUTPARAMS
MDL_CTF_MODEL = xmipp.MDL_CTF_MODEL
MDL_CTF_MODEL2 = xmipp.MDL_CTF_MODEL2
MDL_CTF_SAMPLING_RATE = xmipp.MDL_CTF_SAMPLING_RATE
MDL_CTF_VOLTAGE = xmipp.MDL_CTF_VOLTAGE
MDL_CTF_DEFOCUSU = xmipp.MDL_CTF_DEFOCUSU
MDL_CTF_DEFOCUSV = xmipp.MDL_CTF_DEFOCUSV
MDL_CTF_DEFOCUS_ANGLE = xmipp.MDL_CTF_DEFOCUS_ANGLE
MDL_CTF_DOWNSAMPLE_PERFORMED = xmipp.MDL_CTF_DOWNSAMPLE_PERFORMED
MDL_CTF_CS = xmipp.MDL_CTF_CS
MDL_CTF_CA = xmipp.MDL_CTF_CA
MDL_CTF_GROUP = xmipp.MDL_CTF_GROUP
MDL_CTF_ENERGY_LOSS = xmipp.MDL_CTF_ENERGY_LOSS
MDL_CTF_LENS_STABILITY = xmipp.MDL_CTF_LENS_STABILITY
MDL_CTF_CONVERGENCE_CONE = xmipp.MDL_CTF_CONVERGENCE_CONE
MDL_CTF_LONGITUDINAL_DISPLACEMENT = xmipp.MDL_CTF_LONGITUDINAL_DISPLACEMENT
MDL_CTF_TRANSVERSAL_DISPLACEMENT = xmipp.MDL_CTF_TRANSVERSAL_DISPLACEMENT
MDL_CTF_Q0 = xmipp.MDL_CTF_Q0
MDL_CTF_K = xmipp.MDL_CTF_K
MDL_CTF_BG_GAUSSIAN_K = xmipp.MDL_CTF_BG_GAUSSIAN_K
MDL_CTF_BG_GAUSSIAN_SIGMAU = xmipp.MDL_CTF_BG_GAUSSIAN_SIGMAU
MDL_CTF_BG_GAUSSIAN_SIGMAV = xmipp.MDL_CTF_BG_GAUSSIAN_SIGMAV
MDL_CTF_BG_GAUSSIAN_CU = xmipp.MDL_CTF_BG_GAUSSIAN_CU
MDL_CTF_BG_GAUSSIAN_CV = xmipp.MDL_CTF_BG_GAUSSIAN_CV
MDL_CTF_BG_GAUSSIAN_ANGLE = xmipp.MDL_CTF_BG_GAUSSIAN_ANGLE
MDL_CTF_BG_SQRT_K = xmipp.MDL_CTF_BG_SQRT_K
MDL_CTF_BG_SQRT_U = xmipp.MDL_CTF_BG_SQRT_U
MDL_CTF_BG_SQRT_V = xmipp.MDL_CTF_BG_SQRT_V
MDL_CTF_BG_SQRT_ANGLE = xmipp.MDL_CTF_BG_SQRT_ANGLE
MDL_CTF_BG_BASELINE = xmipp.MDL_CTF_BG_BASELINE
MDL_CTF_BG_GAUSSIAN2_K = xmipp.MDL_CTF_BG_GAUSSIAN2_K
MDL_CTF_BG_GAUSSIAN2_SIGMAU = xmipp.MDL_CTF_BG_GAUSSIAN2_SIGMAU
MDL_CTF_BG_GAUSSIAN2_SIGMAV = xmipp.MDL_CTF_BG_GAUSSIAN2_SIGMAV
MDL_CTF_BG_GAUSSIAN2_CU = xmipp.MDL_CTF_BG_GAUSSIAN2_CU
MDL_CTF_BG_GAUSSIAN2_CV = xmipp.MDL_CTF_BG_GAUSSIAN2_CV
MDL_CTF_BG_GAUSSIAN2_ANGLE = xmipp.MDL_CTF_BG_GAUSSIAN2_ANGLE
MDL_CTF_CRIT_PSDCORRELATION90 = xmipp.MDL_CTF_CRIT_PSDCORRELATION90
MDL_CTF_CRIT_FIRSTZERORATIO = xmipp.MDL_CTF_CRIT_FIRSTZERORATIO
MDL_CTF_CRIT_FIRSTZEROAVG = xmipp.MDL_CTF_CRIT_FIRSTZEROAVG
MDL_CTF_CRIT_FIRSTZERODISAGREEMENT = xmipp.MDL_CTF_CRIT_FIRSTZERODISAGREEMENT
MDL_CTF_CRIT_NORMALITY = xmipp.MDL_CTF_CRIT_NORMALITY
MDL_CTF_CRIT_DAMPING = xmipp.MDL_CTF_CRIT_DAMPING
MDL_CTF_CRIT_PSDRADIALINTEGRAL = xmipp.MDL_CTF_CRIT_PSDRADIALINTEGRAL
MDL_CTF_CRIT_FITTINGSCORE = xmipp.MDL_CTF_CRIT_FITTINGSCORE
MDL_CTF_CRIT_FITTINGCORR13 = xmipp.MDL_CTF_CRIT_FITTINGCORR13
MDL_CTF_CRIT_ICENESS = xmipp.MDL_CTF_CRIT_ICENESS
MDL_CTF_CRIT_PSDVARIANCE = xmipp.MDL_CTF_CRIT_PSDVARIANCE
MDL_CTF_CRIT_PSDPCA1VARIANCE = xmipp.MDL_CTF_CRIT_PSDPCA1VARIANCE
MDL_CTF_CRIT_PSDPCARUNSTEST = xmipp.MDL_CTF_CRIT_PSDPCARUNSTEST
MDL_CTF_PHASE_SHIFT = xmipp.MDL_CTF_PHASE_SHIFT
MDL_CTF_VPP_RADIUS = xmipp.MDL_CTF_VPP_RADIUS
MDL_CUMULATIVE_SSNR = xmipp.MDL_CUMULATIVE_SSNR

MDL_DATATYPE = xmipp.MDL_DATATYPE
MDL_DATE = xmipp.MDL_DATE
MDL_DEFGROUP = xmipp.MDL_DEFGROUP
MDL_DIMENSIONS_3D = xmipp.MDL_DIMENSIONS_3D
MDL_DIMENSIONS_2D = xmipp.MDL_DIMENSIONS_2D
MDL_DM3_IDTAG = xmipp.MDL_DM3_IDTAG
MDL_DM3_NODEID = xmipp.MDL_DM3_NODEID
MDL_DM3_NUMBER_TYPE = xmipp.MDL_DM3_NUMBER_TYPE
MDL_DM3_PARENTID = xmipp.MDL_DM3_PARENTID
MDL_DM3_TAGCLASS = xmipp.MDL_DM3_TAGCLASS
MDL_DM3_TAGNAME = xmipp.MDL_DM3_TAGNAME
MDL_DM3_SIZE = xmipp.MDL_DM3_SIZE
MDL_DM3_VALUE = xmipp.MDL_DM3_VALUE

MDL_ENABLED = xmipp.MDL_ENABLED

MDL_FLIP = xmipp.MDL_FLIP
MDL_FOM = xmipp.MDL_FOM
MDL_FRAME_ID = xmipp.MDL_FRAME_ID

MDL_GATHER_ID = xmipp.MDL_GATHER_ID

MDL_IDX = xmipp.MDL_IDX
MDL_IMAGE = xmipp.MDL_IMAGE
MDL_IMAGE_COVARIANCE = xmipp.MDL_IMAGE_COVARIANCE
MDL_IMAGE_IDX = xmipp.MDL_IMAGE_IDX
MDL_IMAGE_ORIGINAL = xmipp.MDL_IMAGE_ORIGINAL
MDL_IMAGE_REF = xmipp.MDL_IMAGE_REF
MDL_IMAGE_RESIDUAL = xmipp.MDL_IMAGE_RESIDUAL
MDL_IMAGE_TILTED = xmipp.MDL_IMAGE_TILTED
MDL_IMGMD = xmipp.MDL_IMGMD
MDL_IMAGE1 = xmipp.MDL_IMAGE1
MDL_IMAGE2 = xmipp.MDL_IMAGE2
MDL_IMAGE3 = xmipp.MDL_IMAGE3
MDL_IMAGE4 = xmipp.MDL_IMAGE4
MDL_IMAGE5 = xmipp.MDL_IMAGE5
MDL_INTSCALE = xmipp.MDL_INTSCALE
MDL_ITEM_ID = xmipp.MDL_ITEM_ID
MDL_ITER = xmipp.MDL_ITER
MDL_KSTEST = xmipp.MDL_KSTEST
MDL_LL = xmipp.MDL_LL
MDL_MACRO_CMD = xmipp.MDL_MACRO_CMD
MDL_MACRO_CMD_ARGS = xmipp.MDL_MACRO_CMD_ARGS
MDL_MAGNIFICATION = xmipp.MDL_MAGNIFICATION
MDL_MASK = xmipp.MDL_MASK
MDL_MAXCC = xmipp.MDL_MAXCC
MDL_MAX = xmipp.MDL_MAX
MDL_MICROGRAPH = xmipp.MDL_MICROGRAPH
MDL_MICROGRAPH_ID = xmipp.MDL_MICROGRAPH_ID
MDL_MICROGRAPH_MOVIE = xmipp.MDL_MICROGRAPH_MOVIE
MDL_MICROGRAPH_MOVIE_ID = xmipp.MDL_MICROGRAPH_MOVIE_ID
MDL_MICROGRAPH_PARTICLES = xmipp.MDL_MICROGRAPH_PARTICLES
MDL_MICROGRAPH_ORIGINAL = xmipp.MDL_MICROGRAPH_ORIGINAL
MDL_MICROGRAPH_TILTED = xmipp.MDL_MICROGRAPH_TILTED
MDL_MICROGRAPH_TILTED_ORIGINAL = xmipp.MDL_MICROGRAPH_TILTED_ORIGINAL
MDL_MIN = xmipp.MDL_MIN
MDL_MIRRORFRAC = xmipp.MDL_MIRRORFRAC
MDL_MISSINGREGION_NR = xmipp.MDL_MISSINGREGION_NR
MDL_MISSINGREGION_TYPE = xmipp.MDL_MISSINGREGION_TYPE
MDL_MISSINGREGION_THY0 = xmipp.MDL_MISSINGREGION_THY0
MDL_MISSINGREGION_THYF = xmipp.MDL_MISSINGREGION_THYF
MDL_MISSINGREGION_THX0 = xmipp.MDL_MISSINGREGION_THX0
MDL_MISSINGREGION_THXF = xmipp.MDL_MISSINGREGION_THXF
MDL_MLF_CTF = xmipp.MDL_MLF_CTF
MDL_MLF_WIENER = xmipp.MDL_MLF_WIENER
MDL_MLF_SIGNAL = xmipp.MDL_MLF_SIGNAL
MDL_MLF_NOISE = xmipp.MDL_MLF_NOISE
MDL_MODELFRAC = xmipp.MDL_MODELFRAC

MDL_NEIGHBORS = xmipp.MDL_NEIGHBORS
MDL_NEIGHBOR = xmipp.MDL_NEIGHBOR
MDL_NEIGHBORHOOD_RADIUS = xmipp.MDL_NEIGHBORHOOD_RADIUS
MDL_NMA = xmipp.MDL_NMA
MDL_NMA_MODEFILE = xmipp.MDL_NMA_MODEFILE
MDL_NMA_COLLECTIVITY = xmipp.MDL_NMA_COLLECTIVITY
MDL_NMA_MINRANGE = xmipp.MDL_NMA_MINRANGE
MDL_NMA_MAXRANGE = xmipp.MDL_NMA_MAXRANGE
MDL_NMA_SCORE = xmipp.MDL_NMA_SCORE
MDL_NMA_ATOMSHIFT = xmipp.MDL_NMA_ATOMSHIFT
MDL_NOISE_ANGLES = xmipp.MDL_NOISE_ANGLES
MDL_NOISE_PARTICLE_COORD = xmipp.MDL_NOISE_PARTICLE_COORD
MDL_NOISE_COORD = xmipp.MDL_NOISE_COORD
MDL_NOISE_PIXEL_LEVEL = xmipp.MDL_NOISE_PIXEL_LEVEL

MDL_OPTICALFLOW_MEANX = xmipp.MDL_OPTICALFLOW_MEANX
MDL_OPTICALFLOW_MEANY = xmipp.MDL_OPTICALFLOW_MEANY
MDL_OPTICALFLOW_STDX = xmipp.MDL_OPTICALFLOW_STDX
MDL_OPTICALFLOW_STDY = xmipp.MDL_OPTICALFLOW_STDY

MDL_ORDER = xmipp.MDL_ORDER
MDL_ORIGIN_X = xmipp.MDL_ORIGIN_X
MDL_ORIGIN_Y = xmipp.MDL_ORIGIN_Y
MDL_ORIGIN_Z = xmipp.MDL_ORIGIN_Z

MDL_PARTICLE_ID = xmipp.MDL_PARTICLE_ID
MDL_PICKING_AUTOPICKPERCENT = xmipp.MDL_PICKING_AUTOPICKPERCENT
MDL_PICKING_PARTICLE_SIZE = xmipp.MDL_PICKING_PARTICLE_SIZE
MDL_PICKING_STATE = xmipp.MDL_PICKING_STATE
MDL_PICKING_MICROGRAPH_STATE = xmipp.MDL_PICKING_MICROGRAPH_STATE
MDL_PICKING_TEMPLATES = xmipp.MDL_PICKING_TEMPLATES
MDL_PICKING_MANUALPARTICLES_SIZE = xmipp.MDL_PICKING_MANUALPARTICLES_SIZE
MDL_PICKING_AUTOPARTICLES_SIZE = xmipp.MDL_PICKING_AUTOPARTICLES_SIZE
MDL_PMAX = xmipp.MDL_PMAX
MDL_AVGPMAX = xmipp.MDL_AVGPMAX
MDL_PROGRAM = xmipp.MDL_PROGRAM
MDL_PRJ_DIMENSIONS = xmipp.MDL_PRJ_DIMENSIONS
MDL_PRJ_ANGFILE = xmipp.MDL_PRJ_ANGFILE
MDL_PRJ_PSI_NOISE = xmipp.MDL_PRJ_PSI_NOISE
MDL_PRJ_PSI_RANDSTR = xmipp.MDL_PRJ_PSI_RANDSTR
MDL_PRJ_PSI_RANGE = xmipp.MDL_PRJ_PSI_RANGE
MDL_PRJ_ROT_NOISE = xmipp.MDL_PRJ_ROT_NOISE
MDL_PRJ_ROT_RANDSTR = xmipp.MDL_PRJ_ROT_RANDSTR
MDL_PRJ_ROT_RANGE = xmipp.MDL_PRJ_ROT_RANGE
MDL_PRJ_TILT_NOISE = xmipp.MDL_PRJ_TILT_NOISE
MDL_PRJ_TILT_RANDSTR = xmipp.MDL_PRJ_TILT_RANDSTR
MDL_PRJ_TILT_RANGE = xmipp.MDL_PRJ_TILT_RANGE
MDL_PRJ_VOL = xmipp.MDL_PRJ_VOL
MDL_PSD = xmipp.MDL_PSD
MDL_PSD_ENHANCED = xmipp.MDL_PSD_ENHANCED

MDL_RANDOMSEED = xmipp.MDL_RANDOMSEED
MDL_REF3D = xmipp.MDL_REF3D
MDL_REF = xmipp.MDL_REF
MDL_REFMD = xmipp.MDL_REFMD
MDL_RESOLUTION_DPR = xmipp.MDL_RESOLUTION_DPR
MDL_RESOLUTION_ERRORL2 = xmipp.MDL_RESOLUTION_ERRORL2
MDL_RESOLUTION_FRC = xmipp.MDL_RESOLUTION_FRC
MDL_RESOLUTION_FRCRANDOMNOISE = xmipp.MDL_RESOLUTION_FRCRANDOMNOISE
MDL_RESOLUTION_FREQ = xmipp.MDL_RESOLUTION_FREQ
MDL_RESOLUTION_FREQREAL = xmipp.MDL_RESOLUTION_FREQREAL
MDL_RESOLUTION_LOG_STRUCTURE_FACTOR = xmipp.MDL_RESOLUTION_LOG_STRUCTURE_FACTOR
MDL_RESOLUTION_SSNR = xmipp.MDL_RESOLUTION_SSNR
MDL_RESOLUTION_STRUCTURE_FACTOR = xmipp.MDL_RESOLUTION_STRUCTURE_FACTOR
MDL_RESOLUTION_RFACTOR = xmipp.MDL_RESOLUTION_RFACTOR

MDL_SAMPLINGRATE = xmipp.MDL_SAMPLINGRATE
MDL_SAMPLINGRATE_ORIGINAL = xmipp.MDL_SAMPLINGRATE_ORIGINAL
MDL_SAMPLINGRATE_X = xmipp.MDL_SAMPLINGRATE_X
MDL_SAMPLINGRATE_Y = xmipp.MDL_SAMPLINGRATE_Y
MDL_SAMPLINGRATE_Z = xmipp.MDL_SAMPLINGRATE_Z
MDL_SCALE = xmipp.MDL_SCALE
MDL_SELFILE = xmipp.MDL_SELFILE
MDL_SERIE = xmipp.MDL_SERIE
MDL_SHIFT_X = xmipp.MDL_SHIFT_X
MDL_SHIFT_Y = xmipp.MDL_SHIFT_Y
MDL_SHIFT_Z = xmipp.MDL_SHIFT_Z
MDL_SHIFT_X2 = xmipp.MDL_SHIFT_X2
MDL_SHIFT_Y2 = xmipp.MDL_SHIFT_Y2
MDL_SHIFT_X_DIFF = xmipp.MDL_SHIFT_X_DIFF
MDL_SHIFT_Y_DIFF = xmipp.MDL_SHIFT_Y_DIFF
MDL_SHIFT_DIFF = xmipp.MDL_SHIFT_DIFF
MDL_SIGMANOISE = xmipp.MDL_SIGMANOISE
MDL_SIGMAOFFSET = xmipp.MDL_SIGMAOFFSET
MDL_SIGNALCHANGE = xmipp.MDL_SIGNALCHANGE
MDL_STAR_COMMENT = xmipp.MDL_STAR_COMMENT
MDL_STDDEV = xmipp.MDL_STDDEV
MDL_SCORE_BY_ALIGNABILITY_PRECISION = xmipp.MDL_SCORE_BY_ALIGNABILITY_PRECISION
MDL_SCORE_BY_ALIGNABILITY_ACCURACY = xmipp.MDL_SCORE_BY_ALIGNABILITY_ACCURACY
MDL_SCORE_BY_MIRROR = xmipp.MDL_SCORE_BY_MIRROR
MDL_SCORE_BY_ALIGNABILITY_PRECISION_EXP = xmipp.MDL_SCORE_BY_ALIGNABILITY_PRECISION_EXP
MDL_SCORE_BY_ALIGNABILITY_PRECISION_REF = xmipp.MDL_SCORE_BY_ALIGNABILITY_PRECISION_REF
MDL_SCORE_BY_ALIGNABILITY_ACCURACY_EXP = xmipp.MDL_SCORE_BY_ALIGNABILITY_ACCURACY_EXP
MDL_SCORE_BY_ALIGNABILITY_ACCURACY_REF = xmipp.MDL_SCORE_BY_ALIGNABILITY_ACCURACY_REF
MDL_SCORE_BY_ALIGNABILITY_NOISE = xmipp.MDL_SCORE_BY_ALIGNABILITY_NOISE
MDL_SCORE_BY_EMPTINESS = xmipp.MDL_SCORE_BY_EMPTINESS
MDL_SCORE_BY_ENTROPY = xmipp.MDL_SCORE_BY_ENTROPY
MDL_SCORE_BY_GRANULO = xmipp.MDL_SCORE_BY_GRANULO
MDL_SCORE_BY_GINI = xmipp.MDL_SCORE_BY_GINI
MDL_SCORE_BY_LBP = xmipp.MDL_SCORE_BY_LBP
MDL_SCORE_BY_PCA_RESIDUAL=xmipp.MDL_SCORE_BY_PCA_RESIDUAL
MDL_SCORE_BY_PCA_RESIDUAL_PROJ=xmipp.MDL_SCORE_BY_PCA_RESIDUAL_PROJ
MDL_SCORE_BY_PCA_RESIDUAL_EXP=xmipp.MDL_SCORE_BY_PCA_RESIDUAL_EXP
MDL_SCORE_BY_SCREENING = xmipp.MDL_SCORE_BY_SCREENING
MDL_SCORE_BY_VARIANCE = xmipp.MDL_SCORE_BY_VARIANCE
MDL_SCORE_BY_VAR = xmipp.MDL_SCORE_BY_VAR
MDL_SCORE_BY_ZERNIKE = xmipp.MDL_SCORE_BY_ZERNIKE
MDL_SCORE_BY_ZSCORE=xmipp.MDL_SCORE_BY_ZSCORE

MDL_SUM = xmipp.MDL_SUM
MDL_SUMWEIGHT = xmipp.MDL_SUMWEIGHT
MDL_SYMNO = xmipp.MDL_SYMNO

MDL_TIME = xmipp.MDL_TIME
MDL_TRANSFORM_MATRIX = xmipp.MDL_TRANSFORM_MATRIX
MDL_TOMOGRAM_VOLUME = xmipp.MDL_TOMOGRAM_VOLUME
MDL_TOMOGRAMMD = xmipp.MDL_TOMOGRAMMD

MDL_USER = xmipp.MDL_USER

MDL_VOLUME_SCORE_SUM = xmipp.MDL_VOLUME_SCORE_SUM
MDL_VOLUME_SCORE_SUM_TH = xmipp.MDL_VOLUME_SCORE_SUM_TH
MDL_VOLUME_SCORE_MEAN = xmipp.MDL_VOLUME_SCORE_MEAN
MDL_VOLUME_SCORE_MIN = xmipp.MDL_VOLUME_SCORE_MIN
MDL_VOLUME_SCORE1 = xmipp.MDL_VOLUME_SCORE1
MDL_VOLUME_SCORE2 = xmipp.MDL_VOLUME_SCORE2
MDL_VOLUME_SCORE3 = xmipp.MDL_VOLUME_SCORE3
MDL_VOLUME_SCORE4 = xmipp.MDL_VOLUME_SCORE4

MDL_WEIGHT = xmipp.MDL_WEIGHT
MDL_WEIGHT_PRECISION_ALIGNABILITY = xmipp.MDL_WEIGHT_PRECISION_ALIGNABILITY
MDL_WEIGHT_ACCURACY_ALIGNABILITY = xmipp.MDL_WEIGHT_ACCURACY_ALIGNABILITY
MDL_WEIGHT_ALIGNABILITY = xmipp.MDL_WEIGHT_ALIGNABILITY
MDL_WEIGHT_PRECISION_MIRROR = xmipp.MDL_WEIGHT_PRECISION_MIRROR
MDL_WEIGHT_P = xmipp.MDL_WEIGHT_P

MDL_WROBUST = xmipp.MDL_WROBUST

MDL_XCOOR = xmipp.MDL_XCOOR
MDL_XCOOR_TILT = xmipp.MDL_XCOOR_TILT
MDL_XSIZE = xmipp.MDL_XSIZE
MDL_X = xmipp.MDL_X

MDL_YCOOR = xmipp.MDL_YCOOR
MDL_YCOOR_TILT = xmipp.MDL_YCOOR_TILT
MDL_Y = xmipp.MDL_Y
MDL_YSIZE = xmipp.MDL_YSIZE
MDL_ZCOOR = xmipp.MDL_ZCOOR
MDL_Z = xmipp.MDL_Z
MDL_ZSCORE = xmipp.MDL_ZSCORE
MDL_ZSCORE_HISTOGRAM = xmipp.MDL_ZSCORE_HISTOGRAM
MDL_ZSCORE_RESMEAN = xmipp.MDL_ZSCORE_RESMEAN
MDL_ZSCORE_RESVAR = xmipp.MDL_ZSCORE_RESVAR
MDL_ZSCORE_RESCOV = xmipp.MDL_ZSCORE_RESCOV
MDL_ZSCORE_SHAPE1 = xmipp.MDL_ZSCORE_SHAPE1
MDL_ZSCORE_SHAPE2 = xmipp.MDL_ZSCORE_SHAPE2
MDL_ZSCORE_SNR1 = xmipp.MDL_ZSCORE_SNR1
MDL_ZSCORE_SNR2 = xmipp.MDL_ZSCORE_SNR2
MDL_ZSIZE = xmipp.MDL_ZSIZE
MDL_ZSIZE = xmipp.MDL_ZSIZE
MDL_LAST_LABEL = xmipp.MDL_LAST_LABEL

# ----- Label types ------
LABEL_NOTYPE = xmipp.LABEL_NOTYPE
LABEL_INT = xmipp.LABEL_INT
LABEL_BOOL = xmipp.LABEL_BOOL
LABEL_DOUBLE = xmipp.LABEL_DOUBLE
LABEL_VECTOR_DOUBLE = xmipp.LABEL_VECTOR_DOUBLE
LABEL_STRING = xmipp.LABEL_STRING
LABEL_SIZET = xmipp.LABEL_SIZET
LABEL_VECTOR_SIZET = xmipp.LABEL_VECTOR_SIZET
TAGLABEL_NOTAG = xmipp.TAGLABEL_NOTAG
TAGLABEL_TEXTFILE = xmipp.TAGLABEL_TEXTFILE
TAGLABEL_METADATA = xmipp.TAGLABEL_METADATA
TAGLABEL_CTFPARAM = xmipp.TAGLABEL_CTFPARAM
TAGLABEL_IMAGE = xmipp.TAGLABEL_IMAGE
TAGLABEL_VOLUME = xmipp.TAGLABEL_VOLUME
TAGLABEL_STACK = xmipp.TAGLABEL_STACK
TAGLABEL_MICROGRAPH = xmipp.TAGLABEL_MICROGRAPH
TAGLABEL_PSD = xmipp.TAGLABEL_PSD

# ----- RELION labels -------
RLN_AREA_ID = xmipp.RLN_AREA_ID 
RLN_AREA_NAME = xmipp.RLN_AREA_NAME 
RLN_COMMENT = xmipp.RLN_COMMENT 

RLN_CTF_BFACTOR = xmipp.RLN_CTF_BFACTOR 
RLN_CTF_SCALEFACTOR = xmipp.RLN_CTF_SCALEFACTOR 
RLN_CTF_SAMPLING_RATE = xmipp.RLN_CTF_SAMPLING_RATE 
RLN_CTF_VOLTAGE = xmipp.RLN_CTF_VOLTAGE 
RLN_CTF_DEFOCUSU = xmipp.RLN_CTF_DEFOCUSU 
RLN_CTF_DEFOCUSV = xmipp.RLN_CTF_DEFOCUSV 
RLN_CTF_DEFOCUS_ANGLE = xmipp.RLN_CTF_DEFOCUS_ANGLE 
RLN_CTF_CS = xmipp.RLN_CTF_CS 
RLN_CTF_CA = xmipp.RLN_CTF_CA 
RLN_CTF_DETECTOR_PIXEL_SIZE = xmipp.RLN_CTF_DETECTOR_PIXEL_SIZE 
RLN_CTF_ENERGY_LOSS = xmipp.RLN_CTF_ENERGY_LOSS 
RLN_CTF_FOM = xmipp.RLN_CTF_FOM 
RLN_CTF_IMAGE = xmipp.RLN_CTF_IMAGE 
RLN_CTF_LENS_STABILITY = xmipp.RLN_CTF_LENS_STABILITY 
RLN_CTF_MAGNIFICATION = xmipp.RLN_CTF_MAGNIFICATION 
RLN_CTF_CONVERGENCE_CONE = xmipp.RLN_CTF_CONVERGENCE_CONE 
RLN_CTF_LONGITUDINAL_DISPLACEMENT = xmipp.RLN_CTF_LONGITUDINAL_DISPLACEMENT 
RLN_CTF_TRANSVERSAL_DISPLACEMENT = xmipp.RLN_CTF_TRANSVERSAL_DISPLACEMENT 
RLN_CTF_Q0 = xmipp.RLN_CTF_Q0 
RLN_CTF_K = xmipp.RLN_CTF_K 
RLN_CTF_VALUE = xmipp.RLN_CTF_VALUE
RLN_CTF_PHASESHIFT = xmipp.RLN_CTF_PHASESHIFT

RLN_IMAGE_NAME = xmipp.RLN_IMAGE_NAME
RLN_IMAGE_RECONSTRUCT_NAME = xmipp.RLN_IMAGE_RECONSTRUCT_NAME
RLN_IMAGE_ID = xmipp.RLN_IMAGE_ID
RLN_IMAGE_ENABLED = xmipp.RLN_IMAGE_ENABLED
RLN_IMAGE_DATATYPE = xmipp.RLN_IMAGE_DATATYPE
RLN_IMAGE_DIMENSIONALITY = xmipp.RLN_IMAGE_DIMENSIONALITY
RLN_IMAGE_BEAMTILT_X = xmipp.RLN_IMAGE_BEAMTILT_X
RLN_IMAGE_BEAMTILT_Y = xmipp.RLN_IMAGE_BEAMTILT_Y
RLN_IMAGE_BEAMTILT_GROUP = xmipp.RLN_IMAGE_BEAMTILT_GROUP
RLN_IMAGE_COORD_X = xmipp.RLN_IMAGE_COORD_X
RLN_IMAGE_COORD_Y = xmipp.RLN_IMAGE_COORD_Y
RLN_IMAGE_COORD_Z = xmipp.RLN_IMAGE_COORD_Z
RLN_IMAGE_FRAME_NR = xmipp.RLN_IMAGE_FRAME_NR
RLN_IMAGE_MAGNIFICATION_CORRECTION = xmipp.RLN_IMAGE_MAGNIFICATION_CORRECTION
RLN_IMAGE_NORM_CORRECTION = xmipp.RLN_IMAGE_NORM_CORRECTION
RLN_IMAGE_SAMPLINGRATE = xmipp.RLN_IMAGE_SAMPLINGRATE
RLN_IMAGE_SAMPLINGRATE_X = xmipp.RLN_IMAGE_SAMPLINGRATE_X
RLN_IMAGE_SAMPLINGRATE_Y = xmipp.RLN_IMAGE_SAMPLINGRATE_Y
RLN_IMAGE_SAMPLINGRATE_Z = xmipp.RLN_IMAGE_SAMPLINGRATE_Z
RLN_IMAGE_SIZE = xmipp.RLN_IMAGE_SIZE
RLN_IMAGE_SIZEX = xmipp.RLN_IMAGE_SIZEX
RLN_IMAGE_SIZEY = xmipp.RLN_IMAGE_SIZEY
RLN_IMAGE_SIZEZ = xmipp.RLN_IMAGE_SIZEZ
RLN_IMAGE_STATS_MIN = xmipp.RLN_IMAGE_STATS_MIN
RLN_IMAGE_STATS_MAX = xmipp.RLN_IMAGE_STATS_MAX
RLN_IMAGE_STATS_AVG = xmipp.RLN_IMAGE_STATS_AVG
RLN_IMAGE_STATS_STDDEV = xmipp.RLN_IMAGE_STATS_STDDEV
RLN_IMAGE_STATS_SKEW = xmipp.RLN_IMAGE_STATS_SKEW
RLN_IMAGE_STATS_KURT = xmipp.RLN_IMAGE_STATS_KURT
RLN_IMAGE_WEIGHT = xmipp.RLN_IMAGE_WEIGHT

RLN_MATRIX_1_1 = xmipp.RLN_MATRIX_1_1
RLN_MATRIX_1_2 = xmipp.RLN_MATRIX_1_2
RLN_MATRIX_1_3 = xmipp.RLN_MATRIX_1_3
RLN_MATRIX_2_1 = xmipp.RLN_MATRIX_2_1
RLN_MATRIX_2_2 = xmipp.RLN_MATRIX_2_2
RLN_MATRIX_2_3 = xmipp.RLN_MATRIX_2_3
RLN_MATRIX_3_1 = xmipp.RLN_MATRIX_3_1
RLN_MATRIX_3_2 = xmipp.RLN_MATRIX_3_2
RLN_MATRIX_3_3 = xmipp.RLN_MATRIX_3_3

RLN_MICROGRAPH_ID = xmipp.RLN_MICROGRAPH_ID
RLN_MICROGRAPH_MOVIE_NAME = xmipp.RLN_MICROGRAPH_MOVIE_NAME
RLN_MICROGRAPH_NAME = xmipp.RLN_MICROGRAPH_NAME
RLN_MICROGRAPH_TILT_ANGLE = xmipp.RLN_MICROGRAPH_TILT_ANGLE
RLN_MICROGRAPH_TILT_AXIS_DIRECTION = xmipp.RLN_MICROGRAPH_TILT_AXIS_DIRECTION
RLN_MICROGRAPH_TILT_AXIS_OUTOFPLANE = xmipp.RLN_MICROGRAPH_TILT_AXIS_OUTOFPLANE

RLN_MLMODEL_ACCURACY_ROT = xmipp.RLN_MLMODEL_ACCURACY_ROT
RLN_MLMODEL_ACCURACY_TRANS = xmipp.RLN_MLMODEL_ACCURACY_TRANS
RLN_MLMODEL_AVE_PMAX = xmipp.RLN_MLMODEL_AVE_PMAX
RLN_MLMODEL_CURRENT_RESOLUTION = xmipp.RLN_MLMODEL_CURRENT_RESOLUTION
RLN_MLMODEL_CURRENT_SIZE = xmipp.RLN_MLMODEL_CURRENT_SIZE
RLN_MLMODEL_DATA_VS_PRIOR_REF = xmipp.RLN_MLMODEL_DATA_VS_PRIOR_REF
RLN_MLMODEL_DIMENSIONALITY = xmipp.RLN_MLMODEL_DIMENSIONALITY
RLN_MLMODEL_DIMENSIONALITY_DATA = xmipp.RLN_MLMODEL_DIMENSIONALITY_DATA
RLN_MLMODEL_DIFF2_HALVES_REF = xmipp.RLN_MLMODEL_DIFF2_HALVES_REF
RLN_MLMODEL_FSC_HALVES_REF = xmipp.RLN_MLMODEL_FSC_HALVES_REF
RLN_MLMODEL_GROUP_NAME = xmipp.RLN_MLMODEL_GROUP_NAME
RLN_MLMODEL_GROUP_NO = xmipp.RLN_MLMODEL_GROUP_NO
RLN_MLMODEL_GROUP_NR_PARTICLES = xmipp.RLN_MLMODEL_GROUP_NR_PARTICLES
RLN_MLMODEL_GROUP_SCALE_CORRECTION = xmipp.RLN_MLMODEL_GROUP_SCALE_CORRECTION
RLN_MLMODEL_INTERPOLATOR = xmipp.RLN_MLMODEL_INTERPOLATOR
RLN_MLMODEL_LL = xmipp.RLN_MLMODEL_LL
RLN_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION = xmipp.RLN_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION
RLN_MLMODEL_NORM_CORRECTION_AVG = xmipp.RLN_MLMODEL_NORM_CORRECTION_AVG
RLN_MLMODEL_NR_CLASSES = xmipp.RLN_MLMODEL_NR_CLASSES
RLN_MLMODEL_NR_GROUPS = xmipp.RLN_MLMODEL_NR_GROUPS
RLN_MLMODEL_ORIGINAL_SIZE = xmipp.RLN_MLMODEL_ORIGINAL_SIZE
RLN_MLMODEL_ORIENTABILITY_CONTRIBUTION = xmipp.RLN_MLMODEL_ORIENTABILITY_CONTRIBUTION
RLN_MLMODEL_PADDING_FACTOR = xmipp.RLN_MLMODEL_PADDING_FACTOR
RLN_MLMODEL_PDF_CLASS = xmipp.RLN_MLMODEL_PDF_CLASS
RLN_MLMODEL_PRIOR_OFFX_CLASS = xmipp.RLN_MLMODEL_PRIOR_OFFX_CLASS
RLN_MLMODEL_PRIOR_OFFY_CLASS = xmipp.RLN_MLMODEL_PRIOR_OFFY_CLASS
RLN_MLMODEL_PDF_ORIENT = xmipp.RLN_MLMODEL_PDF_ORIENT
RLN_MLMODEL_PIXEL_SIZE = xmipp.RLN_MLMODEL_PIXEL_SIZE
RLN_MLMODEL_POWER_REF = xmipp.RLN_MLMODEL_POWER_REF
RLN_MLMODEL_PRIOR_MODE = xmipp.RLN_MLMODEL_PRIOR_MODE
RLN_MLMODEL_SIGMA_OFFSET = xmipp.RLN_MLMODEL_SIGMA_OFFSET
RLN_MLMODEL_SIGMA_ROT = xmipp.RLN_MLMODEL_SIGMA_ROT
RLN_MLMODEL_SIGMA_TILT = xmipp.RLN_MLMODEL_SIGMA_TILT
RLN_MLMODEL_SIGMA_PSI = xmipp.RLN_MLMODEL_SIGMA_PSI
RLN_MLMODEL_REF_IMAGE = xmipp.RLN_MLMODEL_REF_IMAGE
RLN_MLMODEL_SIGMA2_NOISE = xmipp.RLN_MLMODEL_SIGMA2_NOISE
RLN_MLMODEL_SIGMA2_REF = xmipp.RLN_MLMODEL_SIGMA2_REF
RLN_MLMODEL_SSNR_REF = xmipp.RLN_MLMODEL_SSNR_REF
RLN_MLMODEL_TAU2_FUDGE_FACTOR = xmipp.RLN_MLMODEL_TAU2_FUDGE_FACTOR
RLN_MLMODEL_TAU2_REF = xmipp.RLN_MLMODEL_TAU2_REF
RLN_OPTIMISER_ACCURACY_ROT = xmipp.RLN_OPTIMISER_ACCURACY_ROT
RLN_OPTIMISER_ACCURACY_TRANS = xmipp.RLN_OPTIMISER_ACCURACY_TRANS
RLN_OPTIMISER_ADAPTIVE_FRACTION = xmipp.RLN_OPTIMISER_ADAPTIVE_FRACTION
RLN_OPTIMISER_ADAPTIVE_OVERSAMPLING = xmipp.RLN_OPTIMISER_ADAPTIVE_OVERSAMPLING
RLN_OPTIMISER_AUTO_LOCAL_HP_ORDER = xmipp.RLN_OPTIMISER_AUTO_LOCAL_HP_ORDER
RLN_OPTIMISER_AVAILABLE_MEMORY = xmipp.RLN_OPTIMISER_AVAILABLE_MEMORY
RLN_OPTIMISER_BEST_RESOL_THUS_FAR = xmipp.RLN_OPTIMISER_BEST_RESOL_THUS_FAR
RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS = xmipp.RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS
RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS = xmipp.RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS
RLN_OPTIMISER_CHANGES_OPTIMAL_CLASSES = xmipp.RLN_OPTIMISER_CHANGES_OPTIMAL_CLASSES
RLN_OPTIMISER_COARSE_SIZE = xmipp.RLN_OPTIMISER_COARSE_SIZE
RLN_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED = xmipp.RLN_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED
RLN_OPTIMISER_DATA_STARFILE = xmipp.RLN_OPTIMISER_DATA_STARFILE
RLN_OPTIMISER_DO_AUTO_REFINE = xmipp.RLN_OPTIMISER_DO_AUTO_REFINE
RLN_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES = xmipp.RLN_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES
RLN_OPTIMISER_DO_CORRECT_CTF = xmipp.RLN_OPTIMISER_DO_CORRECT_CTF
RLN_OPTIMISER_DO_CORRECT_MAGNIFICATION = xmipp.RLN_OPTIMISER_DO_CORRECT_MAGNIFICATION
RLN_OPTIMISER_DO_CORRECT_NORM = xmipp.RLN_OPTIMISER_DO_CORRECT_NORM
RLN_OPTIMISER_DO_CORRECT_SCALE = xmipp.RLN_OPTIMISER_DO_CORRECT_SCALE
RLN_OPTIMISER_DO_REALIGN_MOVIES = xmipp.RLN_OPTIMISER_DO_REALIGN_MOVIES
RLN_OPTIMISER_DO_MAP = xmipp.RLN_OPTIMISER_DO_MAP
RLN_OPTIMISER_DO_SOLVENT_FLATTEN = xmipp.RLN_OPTIMISER_DO_SOLVENT_FLATTEN
RLN_OPTIMISER_DO_SKIP_ALIGN = xmipp.RLN_OPTIMISER_DO_SKIP_ALIGN
RLN_OPTIMISER_DO_SKIP_ROTATE = xmipp.RLN_OPTIMISER_DO_SKIP_ROTATE
RLN_OPTIMISER_DO_SPLIT_RANDOM_HALVES = xmipp.RLN_OPTIMISER_DO_SPLIT_RANDOM_HALVES
RLN_OPTIMISER_DO_ZERO_MASK = xmipp.RLN_OPTIMISER_DO_ZERO_MASK
RLN_OPTIMISER_FIX_SIGMA_NOISE = xmipp.RLN_OPTIMISER_FIX_SIGMA_NOISE
RLN_OPTIMISER_FIX_SIGMA_OFFSET = xmipp.RLN_OPTIMISER_FIX_SIGMA_OFFSET
RLN_OPTIMISER_FIX_TAU = xmipp.RLN_OPTIMISER_FIX_TAU
RLN_OPTIMISER_HAS_CONVERGED = xmipp.RLN_OPTIMISER_HAS_CONVERGED
RLN_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT = xmipp.RLN_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT
RLN_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO = xmipp.RLN_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO
RLN_OPTIMISER_HIGHRES_LIMIT_EXP = xmipp.RLN_OPTIMISER_HIGHRES_LIMIT_EXP
RLN_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK = xmipp.RLN_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK
RLN_OPTIMISER_INCR_SIZE = xmipp.RLN_OPTIMISER_INCR_SIZE
RLN_OPTIMISER_ITERATION_NO = xmipp.RLN_OPTIMISER_ITERATION_NO
RLN_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES = xmipp.RLN_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES
RLN_OPTIMISER_MAGNIFICATION_RANGE = xmipp.RLN_OPTIMISER_MAGNIFICATION_RANGE
RLN_OPTIMISER_MAGNIFICATION_STEP = xmipp.RLN_OPTIMISER_MAGNIFICATION_STEP
RLN_OPTIMISER_MAX_COARSE_SIZE = xmipp.RLN_OPTIMISER_MAX_COARSE_SIZE
RLN_OPTIMISER_MAX_NR_POOL = xmipp.RLN_OPTIMISER_MAX_NR_POOL
RLN_OPTIMISER_MODEL_STARFILE = xmipp.RLN_OPTIMISER_MODEL_STARFILE
RLN_OPTIMISER_MODEL_STARFILE2 = xmipp.RLN_OPTIMISER_MODEL_STARFILE2
RLN_OPTIMISER_NR_ITERATIONS = xmipp.RLN_OPTIMISER_NR_ITERATIONS
RLN_OPTIMISER_NR_ITER_WO_RESOL_GAIN = xmipp.RLN_OPTIMISER_NR_ITER_WO_RESOL_GAIN
RLN_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES = xmipp.RLN_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES
RLN_OPTIMISER_OUTPUT_ROOTNAME = xmipp.RLN_OPTIMISER_OUTPUT_ROOTNAME
RLN_OPTIMISER_PARTICLE_DIAMETER = xmipp.RLN_OPTIMISER_PARTICLE_DIAMETER
RLN_OPTIMISER_RADIUS_MASK_3D_MAP = xmipp.RLN_OPTIMISER_RADIUS_MASK_3D_MAP
RLN_OPTIMISER_RADIUS_MASK_EXP_PARTICLES = xmipp.RLN_OPTIMISER_RADIUS_MASK_EXP_PARTICLES
RLN_OPTIMISER_RANDOM_SEED = xmipp.RLN_OPTIMISER_RANDOM_SEED
RLN_OPTIMISER_REFS_ARE_CTF_CORRECTED = xmipp.RLN_OPTIMISER_REFS_ARE_CTF_CORRECTED
RLN_OPTIMISER_SAMPLING_STARFILE = xmipp.RLN_OPTIMISER_SAMPLING_STARFILE
RLN_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES = xmipp.RLN_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES
RLN_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS = xmipp.RLN_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS
RLN_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS = xmipp.RLN_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS
RLN_OPTIMISER_SOLVENT_MASK_NAME = xmipp.RLN_OPTIMISER_SOLVENT_MASK_NAME
RLN_OPTIMISER_SOLVENT_MASK2_NAME = xmipp.RLN_OPTIMISER_SOLVENT_MASK2_NAME
RLN_OPTIMISER_TAU_SPECTRUM_NAME = xmipp.RLN_OPTIMISER_TAU_SPECTRUM_NAME
RLN_OPTIMISER_USE_TOO_COARSE_SAMPLING = xmipp.RLN_OPTIMISER_USE_TOO_COARSE_SAMPLING
RLN_OPTIMISER_WIDTH_MASK_EDGE = xmipp.RLN_OPTIMISER_WIDTH_MASK_EDGE
RLN_ORIENT_FLIP = xmipp.RLN_ORIENT_FLIP
RLN_ORIENT_ID = xmipp.RLN_ORIENT_ID
RLN_ORIENT_ORIGIN_X = xmipp.RLN_ORIENT_ORIGIN_X
RLN_ORIENT_ORIGIN_X_PRIOR = xmipp.RLN_ORIENT_ORIGIN_X_PRIOR
RLN_ORIENT_ORIGIN_Y = xmipp.RLN_ORIENT_ORIGIN_Y
RLN_ORIENT_ORIGIN_Y_PRIOR = xmipp.RLN_ORIENT_ORIGIN_Y_PRIOR
RLN_ORIENT_ORIGIN_Z = xmipp.RLN_ORIENT_ORIGIN_Z
RLN_ORIENT_ORIGIN_Z_PRIOR = xmipp.RLN_ORIENT_ORIGIN_Z_PRIOR
RLN_ORIENT_ROT = xmipp.RLN_ORIENT_ROT
RLN_ORIENT_ROT_PRIOR = xmipp.RLN_ORIENT_ROT_PRIOR
RLN_ORIENT_TILT = xmipp.RLN_ORIENT_TILT
RLN_ORIENT_TILT_PRIOR = xmipp.RLN_ORIENT_TILT_PRIOR
RLN_ORIENT_PSI = xmipp.RLN_ORIENT_PSI
RLN_ORIENT_PSI_PRIOR = xmipp.RLN_ORIENT_PSI_PRIOR
RLN_ORIENT_PSI_PRIOR_FLIP_RATIO = xmipp.RLN_ORIENT_PSI_PRIOR_FLIP_RATIO
RLN_PARTICLE_AUTOPICK_FOM = xmipp.RLN_PARTICLE_AUTOPICK_FOM
RLN_PARTICLE_CLASS = xmipp.RLN_PARTICLE_CLASS
RLN_PARTICLE_DLL = xmipp.RLN_PARTICLE_DLL
RLN_PARTICLE_ID = xmipp.RLN_PARTICLE_ID
RLN_PARTICLE_FOM = xmipp.RLN_PARTICLE_FOM
RLN_PARTICLE_KL_DIVERGENCE = xmipp.RLN_PARTICLE_KL_DIVERGENCE
RLN_PARTICLE_MOVIE_RUNNING_AVG = xmipp.RLN_PARTICLE_MOVIE_RUNNING_AVG
RLN_PARTICLE_RANDOM_SUBSET = xmipp.RLN_PARTICLE_RANDOM_SUBSET
RLN_PARTICLE_NAME = xmipp.RLN_PARTICLE_NAME
RLN_PARTICLE_ORI_NAME = xmipp.RLN_PARTICLE_ORI_NAME
RLN_PARTICLE_NR_SIGNIFICANT_SAMPLES = xmipp.RLN_PARTICLE_NR_SIGNIFICANT_SAMPLES
RLN_PARTICLE_NR_FRAMES = xmipp.RLN_PARTICLE_NR_FRAMES
RLN_PARTICLE_NR_FRAMES_AVG = xmipp.RLN_PARTICLE_NR_FRAMES_AVG
RLN_PARTICLE_PMAX = xmipp.RLN_PARTICLE_PMAX

# New helical labes in Relion 2.x
RLN_MLMODEL_HELICAL_NR_ASU = xmipp.RLN_MLMODEL_HELICAL_NR_ASU
RLN_MLMODEL_HELICAL_TWIST = xmipp.RLN_MLMODEL_HELICAL_TWIST
RLN_MLMODEL_HELICAL_TWIST_MIN = xmipp.RLN_MLMODEL_HELICAL_TWIST_MIN	
RLN_MLMODEL_HELICAL_TWIST_MAX = xmipp.RLN_MLMODEL_HELICAL_TWIST_MAX
RLN_MLMODEL_HELICAL_TWIST_INITIAL_STEP = xmipp.RLN_MLMODEL_HELICAL_TWIST_INITIAL_STEP
RLN_MLMODEL_HELICAL_RISE = xmipp.RLN_MLMODEL_HELICAL_RISE
RLN_MLMODEL_HELICAL_RISE_MIN = xmipp.RLN_MLMODEL_HELICAL_RISE_MIN
RLN_MLMODEL_HELICAL_RISE_MAX = xmipp.RLN_MLMODEL_HELICAL_RISE_MAX
RLN_MLMODEL_HELICAL_RISE_INITIAL_STEP = xmipp.RLN_MLMODEL_HELICAL_RISE_INITIAL_STEP
RLN_OPTIMISER_DO_HELICAL_REFINE = xmipp.RLN_OPTIMISER_DO_HELICAL_REFINE
RLN_OPTIMISER_HELICAL_TWIST_INITIAL = xmipp.RLN_OPTIMISER_HELICAL_TWIST_INITIAL
RLN_OPTIMISER_HELICAL_RISE_INITIAL = xmipp.RLN_OPTIMISER_HELICAL_RISE_INITIAL
RLN_OPTIMISER_HELICAL_Z_PERCENTAGE = xmipp.RLN_OPTIMISER_HELICAL_Z_PERCENTAGE
RLN_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER = xmipp.RLN_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER
RLN_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER = xmipp.RLN_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER
RLN_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT = xmipp.RLN_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT
RLN_OPTIMISER_HELICAL_SIGMA_DISTANCE = xmipp.RLN_OPTIMISER_HELICAL_SIGMA_DISTANCE
RLN_OPTIMISER_IGNORE_HELICAL_SYMMETRY = xmipp.RLN_OPTIMISER_IGNORE_HELICAL_SYMMETRY
RLN_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED = xmipp.RLN_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED
RLN_PARTICLE_HELICAL_TUBE_ID = xmipp.RLN_PARTICLE_HELICAL_TUBE_ID
RLN_PARTICLE_HELICAL_TUBE_PITCH = xmipp.RLN_PARTICLE_HELICAL_TUBE_PITCH
RLN_PARTICLE_HELICAL_TRACK_LENGTH = xmipp.RLN_PARTICLE_HELICAL_TRACK_LENGTH
RLN_SAMPLING_HELICAL_OFFSET_STEP = xmipp.RLN_SAMPLING_HELICAL_OFFSET_STEP

RLN_POSTPROCESS_BFACTOR = xmipp.RLN_POSTPROCESS_BFACTOR
RLN_POSTPROCESS_FINAL_RESOLUTION = xmipp.RLN_POSTPROCESS_FINAL_RESOLUTION
RLN_POSTPROCESS_FSC_TRUE = xmipp.RLN_POSTPROCESS_FSC_TRUE
RLN_POSTPROCESS_FSC_MASKED = xmipp.RLN_POSTPROCESS_FSC_MASKED
RLN_POSTPROCESS_FSC_UNMASKED = xmipp.RLN_POSTPROCESS_FSC_UNMASKED
RLN_POSTPROCESS_FSC_RANDOM_MASKED = xmipp.RLN_POSTPROCESS_FSC_RANDOM_MASKED
RLN_POSTPROCESS_GUINIER_FIT_CORRELATION = xmipp.RLN_POSTPROCESS_GUINIER_FIT_CORRELATION
RLN_POSTPROCESS_GUINIER_FIT_INTERCEPT = xmipp.RLN_POSTPROCESS_GUINIER_FIT_INTERCEPT
RLN_POSTPROCESS_GUINIER_FIT_SLOPE = xmipp.RLN_POSTPROCESS_GUINIER_FIT_SLOPE
RLN_POSTPROCESS_GUINIER_VALUE_IN = xmipp.RLN_POSTPROCESS_GUINIER_VALUE_IN
RLN_POSTPROCESS_GUINIER_VALUE_INVMTF = xmipp.RLN_POSTPROCESS_GUINIER_VALUE_INVMTF
RLN_POSTPROCESS_GUINIER_VALUE_WEIGHTED = xmipp.RLN_POSTPROCESS_GUINIER_VALUE_WEIGHTED
RLN_POSTPROCESS_GUINIER_VALUE_SHARPENED = xmipp.RLN_POSTPROCESS_GUINIER_VALUE_SHARPENED
RLN_POSTPROCESS_GUINIER_VALUE_INTERCEPT = xmipp.RLN_POSTPROCESS_GUINIER_VALUE_INTERCEPT
RLN_POSTPROCESS_GUINIER_RESOL_SQUARED = xmipp.RLN_POSTPROCESS_GUINIER_RESOL_SQUARED
RLN_POSTPROCESS_MTF_VALUE = xmipp.RLN_POSTPROCESS_MTF_VALUE 
RLN_SAMPLING_IS_3D = xmipp.RLN_SAMPLING_IS_3D
RLN_SAMPLING_IS_3D_TRANS = xmipp.RLN_SAMPLING_IS_3D_TRANS
RLN_SAMPLING_HEALPIX_ORDER = xmipp.RLN_SAMPLING_HEALPIX_ORDER
RLN_SAMPLING_LIMIT_TILT = xmipp.RLN_SAMPLING_LIMIT_TILT
RLN_SAMPLING_OFFSET_RANGE = xmipp.RLN_SAMPLING_OFFSET_RANGE
RLN_SAMPLING_OFFSET_STEP = xmipp.RLN_SAMPLING_OFFSET_STEP
RLN_SAMPLING_PERTURB = xmipp.RLN_SAMPLING_PERTURB
RLN_SAMPLING_PERTURBATION_FACTOR = xmipp.RLN_SAMPLING_PERTURBATION_FACTOR
RLN_SAMPLING_PRIOR_MODE = xmipp.RLN_SAMPLING_PRIOR_MODE
RLN_SAMPLING_PSI_STEP = xmipp.RLN_SAMPLING_PSI_STEP
RLN_SAMPLING_SIGMA_ROT = xmipp.RLN_SAMPLING_SIGMA_ROT
RLN_SAMPLING_SIGMA_TILT = xmipp.RLN_SAMPLING_SIGMA_TILT
RLN_SAMPLING_SIGMA_PSI = xmipp.RLN_SAMPLING_SIGMA_PSI
RLN_SAMPLING_SYMMETRY = xmipp.RLN_SAMPLING_SYMMETRY

RLN_SELECTED = xmipp.RLN_SELECTED
RLN_SELECT_PARTICLES_ZSCORE = xmipp.RLN_SELECT_PARTICLES_ZSCORE
RLN_SORTED_IDX = xmipp.RLN_SORTED_IDX
RLN_PERFRAME_CUMULATIVE_WEIGHT = xmipp.RLN_PERFRAME_CUMULATIVE_WEIGHT
RLN_PERFRAME_RELATIVE_WEIGHT = xmipp.RLN_PERFRAME_RELATIVE_WEIGHT

RLN_RESOLUTION = xmipp.RLN_RESOLUTION
RLN_RESOLUTION_ANGSTROM = xmipp.RLN_RESOLUTION_ANGSTROM
RLN_RESOLUTION_INVPIXEL = xmipp.RLN_RESOLUTION_INVPIXEL
RLN_SPECTRAL_IDX = xmipp.RLN_SPECTRAL_IDX

# new labels in Relion 2.1
RLN_MLMODEL_ESTIM_RESOL_REF = xmipp.RLN_MLMODEL_ESTIM_RESOL_REF
RLN_MLMODEL_FOURIER_COVERAGE_REF = xmipp.RLN_MLMODEL_FOURIER_COVERAGE_REF
RLN_MLMODEL_FOURIER_COVERAGE_TOTAL_REF = xmipp.RLN_MLMODEL_FOURIER_COVERAGE_TOTAL_REF
RLN_OPTIMISER_LOCAL_SYMMETRY_FILENAME = xmipp.RLN_OPTIMISER_LOCAL_SYMMETRY_FILENAME

# SGD labels in Relion 2.1
RLN_MLMODEL_SGD_GRADIENT_IMAGE = xmipp.RLN_MLMODEL_SGD_GRADIENT_IMAGE
RLN_OPTIMISER_DO_SGD = xmipp.RLN_OPTIMISER_DO_SGD
RLN_OPTIMISER_SGD_MU = xmipp.RLN_OPTIMISER_SGD_MU
RLN_OPTIMISER_SGD_SIGMA2FUDGE_INI = xmipp.RLN_OPTIMISER_SGD_SIGMA2FUDGE_INI
RLN_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE = xmipp.RLN_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE
RLN_OPTIMISER_SGD_SUBSET_START = xmipp.RLN_OPTIMISER_SGD_SUBSET_START
RLN_OPTIMISER_SGD_SUBSET_SIZE = xmipp.RLN_OPTIMISER_SGD_SUBSET_SIZE
RLN_OPTIMISER_SGD_WRITE_EVERY_SUBSET = xmipp.RLN_OPTIMISER_SGD_WRITE_EVERY_SUBSET
RLN_OPTIMISER_SGD_MAX_SUBSETS = xmipp.RLN_OPTIMISER_SGD_MAX_SUBSETS
RLN_OPTIMISER_SGD_STEPSIZE = xmipp.RLN_OPTIMISER_SGD_STEPSIZE
RLN_OPTIMISER_HIGHRES_LIMIT_SGD = xmipp.RLN_OPTIMISER_HIGHRES_LIMIT_SGD

# dataTypes constants
DT_DEFAULT = xmipp.DT_DEFAULT
DT_UNKNOWN = xmipp.DT_UNKNOWN
DT_UCHAR = xmipp.DT_UCHAR
DT_SCHAR = xmipp.DT_SCHAR
DT_USHORT = xmipp.DT_USHORT
DT_SHORT = xmipp.DT_SHORT
DT_UINT = xmipp.DT_UINT
DT_INT = xmipp.DT_INT
DT_LONG = xmipp.DT_LONG
DT_FLOAT = xmipp.DT_FLOAT
DT_DOUBLE = xmipp.DT_DOUBLE
DT_COMPLEXSHORT = xmipp.DT_COMPLEXSHORT
DT_COMPLEXINT = xmipp.DT_COMPLEXINT
DT_COMPLEXFLOAT = xmipp.DT_COMPLEXFLOAT
DT_COMPLEXDOUBLE = xmipp.DT_COMPLEXDOUBLE

LABEL_TYPES = {
               LABEL_SIZET: long,
               LABEL_DOUBLE: float,
               LABEL_INT: int,
               LABEL_BOOL: bool              
               }

