# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se)
# *
# * SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from collections import OrderedDict

import pwem.emlib.metadata as md

# ----------------- Constants values ------------------------------------------

RELION_HOME = 'RELION_HOME'
RELION_CUDA_LIB = 'RELION_CUDA_LIB'

# Supported versions:
V3_0 = '3.0'
V3_1 = '3.1'

MASK_FILL_ZERO = 0
MASK_FILL_NOISE = 1

ANGULAR_SAMPLING_LIST = ['30', '15', '7.5', '3.7', '1.8', '0.9', '0.5', '0.2', '0.1']

# Micrograph type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1

# Picking protocols
REF_AVERAGES = 0
REF_VOLUME = REF_BLOBS = 1

RUN_OPTIMIZE = 0  # Run only on several micrographs to optimize parameters
RUN_COMPUTE = 1  # Run the picking for all micrographs after optimize

MICS_AUTO = 0
MICS_SUBSET = 1

# Used in ctf-refinment
FIT_NO = 0
FIT_MICS = 1
FIT_PARTS = 2


# Protocol create mask 3d
MASK_AND = 0
MASK_OR = 1
MASK_AND_NOT = 2
MASK_OR_NOT = 3

# Protocol export particles
STACK_NONE = 0
STACK_MULT = 1
STACK_ONE = 2

# Viewer constants
ITER_LAST = 0
ITER_SELECTION = 1

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

CHIMERADATAVIEW = 0

CLASSES_ALL = 0
CLASSES_SEL = 1

FSC_CORRECTED = 0
FSC_UNMASKEDMAPS = 1
FSC_MASKEDMAPS = 2
FSC_RANDOMIZED = 3
FSC_ALL = 4

# Color maps
COLOR_JET = 0
COLOR_TERRAIN = 1
COLOR_GIST_EARTH = 2
COLOR_GIST_NCAR = 3
COLOR_GNU_PLOT = 4
COLOR_GNU_PLOT2 = 5
COLOR_OTHER = 6

COLOR_CHOICES = OrderedDict()

COLOR_CHOICES[COLOR_JET] = 'jet'
COLOR_CHOICES[COLOR_TERRAIN] = 'terrain'
COLOR_CHOICES[COLOR_GIST_EARTH] = 'gist_earth'
COLOR_CHOICES[COLOR_GIST_NCAR] = 'gist_ncar'
COLOR_CHOICES[COLOR_GNU_PLOT] = 'gnuplot'
COLOR_CHOICES[COLOR_GNU_PLOT2] = 'gnuplot2'
COLOR_CHOICES[COLOR_OTHER] = 'other'

# Axis code
AX_X = 0
AX_Y = 1
AX_Z = 2

# Constant to define the behaviour when
# dealing with binaries conversion

# Only convert if the format is not supported
# create link otherwise
CONVERT_IF_NEEDED = 0
CONVERT_ALWAYS = 1  # Always convert
CONVERT_FILENAME = 2  # Only convert the filename

# Map from Xmipp labels to Relion labels names
XMIPP_RELION_LABELS = {
    md.MDL_ANGLE_ROT: 'rlnAngleRot',
    md.MDL_ANGLE_TILT: 'rlnAngleTilt',
    md.MDL_ANGLE_PSI: 'rlnAnglePsi',
    md.MDL_AVG_CHANGES_ORIENTATIONS: 'rlnChangesOptimalOrientations',
    md.MDL_AVG_CHANGES_OFFSETS: 'rlnChangesOptimalOffsets',
    md.MDL_AVG_CHANGES_CLASSES: 'rlnChangesOptimalClasses',
    md.MDL_AVGPMAX: 'rlnAveragePmax',
    md.MDL_CLASS_PERCENTAGE: 'rlnClassDistribution',
    md.MDL_CTF_CA: 'rlnChromaticAberration',
    md.MDL_CTF_CONVERGENCE_CONE: 'rlnConvergenceCone',
    md.MDL_CTF_CS: 'rlnSphericalAberration',
    md.MDL_CTF_DEFOCUSU: 'rlnDefocusU',
    md.MDL_CTF_DEFOCUSV: 'rlnDefocusV',
    md.MDL_CTF_DEFOCUS_ANGLE: 'rlnDefocusAngle',
    md.MDL_CTF_ENERGY_LOSS: 'rlnEnergyLoss',
    md.MDL_CTF_LENS_STABILITY: 'rlnLensStability',
    md.MDL_CTF_LONGITUDINAL_DISPLACEMENT: 'rlnLongitudinalDisplacement',
    md.MDL_CTF_TRANSVERSAL_DISPLACEMENT: 'rlnTransversalDisplacement',
    md.MDL_CTF_Q0: 'rlnAmplitudeContrast',
    md.MDL_CTF_SAMPLING_RATE: 'rlnDetectorPixelSize',
    md.MDL_CTF_VOLTAGE: 'rlnVoltage',
    md.MDL_DATATYPE: 'rlnDataType',
    md.MDL_DEFGROUP: 'rlnGroupNumber',
    md.MDL_ENABLED: 'rlnEnabled',
    md.MDL_IMAGE: 'rlnImageName',
    md.MDL_IMAGE_REF: 'rlnReferenceImage',
    md.MDL_ITEM_ID: 'rlnImageId',
    md.MDL_MAXCC: 'rlnLogLikelihood',
    md.MDL_PSD: 'rlnCtfImage',  # relion 1.3
    md.MDL_LL: 'rlnLogLikeliContribution',
    md.MDL_MAGNIFICATION: 'rlnMagnification',
    md.MDL_MICROGRAPH: 'rlnMicrographName',
    md.MDL_MICROGRAPH_ID: 'rlnMicrographId',
    md.MDL_REF: 'rlnClassNumber',
    md.MDL_RESOLUTION_FREQREAL: 'rlnAngstromResolution',
    md.MDL_RESOLUTION_FRC: 'rlnGoldStandardFsc',
    md.MDL_RESOLUTION_FREQ: 'rlnResolution',
    md.MDL_RESOLUTION_SSNR: 'rlnSsnrMap',
    md.MDL_RANDOMSEED: 'rlnRandomSubset',
    md.MDL_SAMPLINGRATE: 'rlnPixelSize',
    md.MDL_SAMPLINGRATE_ORIGINAL: 'rlnPixelSize',
    md.MDL_SCALE: 'rlnMagnificationCorrection',
    md.MDL_SHIFT_X: 'rlnOriginX',
    md.MDL_SHIFT_Y: 'rlnOriginY',
    md.MDL_SHIFT_Z: 'rlnOriginZ',
    md.MDL_PMAX: 'rlnMaxValueProbDistribution',
    md.MDL_SAMPLINGRATE_X: 'rlnSamplingRateX',
    md.MDL_SAMPLINGRATE_Y: 'rlnSamplingRateY',
    md.MDL_SAMPLINGRATE_Z: 'rlnSamplingRateZ',
    md.MDL_XCOOR: 'rlnCoordinateX',
    md.MDL_XSIZE: 'rlnImageSizeX',
    md.MDL_YCOOR: 'rlnCoordinateY',
    md.MDL_YSIZE: 'rlnImageSizeY',
    md.MDL_WEIGHT: 'rlnNrOfSignificantSamples',
    md.MDL_ZSIZE: 'rlnImageSizeZ',
    # relion 1.3
    md.MDL_IMAGE2: 'rlnParticleName',
    md.MDL_IMAGE_ORIGINAL: 'rlnOriginalParticleName',
    md.MDL_SERIE: 'rlnGroupName'
}

XMIPP_RELION_LABELS_EXTRA = {
    # Following labels have no correct matching, just to
    # pick one with the same datatype
    md.MDL_ANGLE_Y2: 'rlnOrientationalPriorMode',  # int
    md.MDL_BLOCK_NUMBER: 'rlnGroupNumber',  # just one
    md.MDL_COUNT: 'rlnGroupNrParticles',  # just one
    md.MDL_CTF_CRIT_FITTINGSCORE: 'rlnCtfFigureOfMerit',  # just one
    md.MDL_CTF_CRIT_NORMALITY: 'rlnNormCorrection',  # just one
    md.MDL_CTF_CRIT_PSDVARIANCE: 'rlnCtfValue',  # just one
    md.MDL_CTF_CRIT_PSDCORRELATION90: 'rlnCtfBfactor',  # just one
    md.MDL_CRYSTAL_CELLX: 'rlnReferenceDimensionality',
    md.MDL_CRYSTAL_CELLY: 'rlnOriginalImageSize',
    md.MDL_CRYSTAL_DISAPPEAR_THRE: 'rlnCurrentResolution',
    md.MDL_PICKING_PARTICLE_SIZE: 'rlnCurrentImageSize',  # int
    md.MDL_PICKING_AUTOPICKPERCENT: 'rlnPaddingFactor',  # int
    md.MDL_PICKING_TEMPLATES: 'rlnFourierSpaceInterpolator',  # int
    md.MDL_COLOR: 'rlnMinRadiusNnInterpolation',  # int
    md.MDL_DM3_NUMBER_TYPE: 'rlnNrClasses',  # int
    md.MDL_DM3_SIZE: 'rlnNrGroups',  # int
    md.MDL_NMA_COLLECTIVITY: 'rlnTau2FudgeFactor',  # double
    md.MDL_NMA_SCORE: 'rlnNormCorrectionAverage',  # double
    md.MDL_SIGMAOFFSET: 'rlnSigmaOffsets',  # double
    md.MDL_MLF_WIENER: 'rlnOrientationDistribution',  # double
    md.MDL_IDX: 'rlnSpectralIndex',  # int
    md.MDL_MLF_NOISE: 'rlnSigma2Noise',  # double
    md.MDL_DM3_TAGNAME: 'rlnGroupName',  # string
    md.MDL_MLF_SIGNAL: 'rlnGroupScaleCorrection',  # double

    md.MDL_FOM: 'rlnAutopickFigureOfMerit',
    md.MDL_ZSCORE_SHAPE1: 'rlnAccuracyRotations',
    md.MDL_ZSCORE_SHAPE2: 'rlnAccuracyTranslations',
    md.MDL_ZSCORE_SNR1: 'rlnReferenceSigma2',
    md.MDL_ZSCORE_SNR2: 'rlnReferenceTau2',
    md.MDL_ZSCORE_RESCOV: 'rlnSpectralOrientabilityContribution',

    # Keep relion 1.3 new labels at the end
    # just in case we want to keep 1.2 compatibility
    md.MDL_ANGLE_ROT2: 'rlnAngleRotPrior',
    md.MDL_ANGLE_TILT2: 'rlnAngleTiltPrior',
    md.MDL_ANGLE_PSI2: 'rlnAnglePsiPrior',
    md.MDL_SHIFT_X2: 'rlnOriginXPrior',
    md.MDL_SHIFT_Y2: 'rlnOriginYPrior',
    md.MDL_ZSCORE: 'rlnParticleSelectZScore',
    md.MDL_COUNT2: 'rlnNrOfFrames',
    # Not the best labels, but just to grab some
    md.MDL_CLASSIFICATION_DPR_05: 'rlnClassPriorOffsetX',
    md.MDL_CLASSIFICATION_FRC_05: 'rlnClassPriorOffsetY',
    # md.MDL_CLASSIFICATION_INTRACLASS_DISTANCE: ''
}

# This dictionary will be used to map
# between CTFModel properties and Xmipp labels

ACQUISITION_DICT = OrderedDict([
    ("_amplitudeContrast", md.RLN_CTF_Q0),
    ("_sphericalAberration", md.RLN_CTF_CS),
    ("_voltage", md.RLN_CTF_VOLTAGE),
    ("_magnification", md.RLN_CTF_MAGNIFICATION)
])

COOR_DICT = OrderedDict([
    ("_x", md.RLN_IMAGE_COORD_X),
    ("_y", md.RLN_IMAGE_COORD_Y)
])

COOR_EXTRA_LABELS = [
    # Additional autopicking-related metadata
    md.RLN_PARTICLE_AUTOPICK_FOM,
    md.RLN_PARTICLE_CLASS,
    md.RLN_ORIENT_PSI
]

CTF_DICT = OrderedDict([
    ("_defocusU", md.RLN_CTF_DEFOCUSU),
    ("_defocusV", md.RLN_CTF_DEFOCUSV),
    ("_defocusAngle", md.RLN_CTF_DEFOCUS_ANGLE)
])

CTF_PSD_DICT = OrderedDict([
    ("_psdFile", md.RLN_CTF_IMAGE)
])

CTF_EXTRA_LABELS = [
    md.RLN_CTF_FOM,
    md.RLN_CTF_PHASESHIFT,
    # In Relion the ctf also contains acquisition information
    md.RLN_CTF_Q0,
    md.RLN_CTF_CS,
    md.RLN_CTF_VOLTAGE,
    md.RLN_CTF_MAGNIFICATION,
    md.RLN_CTF_DETECTOR_PIXEL_SIZE
]

# Some extra labels
IMAGE_EXTRA_LABELS = [
    md.RLN_SELECT_PARTICLES_ZSCORE,
    md.RLN_IMAGE_FRAME_NR
]

# Extra labels for movie refinement & polishing
MOVIE_EXTRA_LABELS = [
    md.RLN_PARTICLE_NR_FRAMES,
    md.RLN_PARTICLE_NR_FRAMES_AVG,
    md.RLN_PARTICLE_MOVIE_RUNNING_AVG,
    md.RLN_PARTICLE_ORI_NAME,
    md.RLN_MLMODEL_GROUP_NAME,
    # the following is required for polishing to work
    md.RLN_PARTICLE_DLL,
    md.RLN_PARTICLE_PMAX,
    md.RLN_IMAGE_NORM_CORRECTION,
    md.RLN_PARTICLE_NR_SIGNIFICANT_SAMPLES,
    md.RLN_PARTICLE_RANDOM_SUBSET,
    md.RLN_ORIENT_ORIGIN_X_PRIOR,
    md.RLN_ORIENT_ORIGIN_Y_PRIOR,
    md.RLN_ORIENT_PSI_PRIOR,
    md.RLN_ORIENT_ROT_PRIOR,
    md.RLN_ORIENT_TILT_PRIOR
]

ALIGNMENT_DICT = OrderedDict([
    ("_rlnOriginX", md.RLN_ORIENT_ORIGIN_X),
    ("_rlnOriginY", md.RLN_ORIENT_ORIGIN_Y),
    ("_rlnOriginZ", md.RLN_ORIENT_ORIGIN_Z),
    ("_rlnAngleRot", md.RLN_ORIENT_ROT),
    ("_rlnAngleTilt", md.RLN_ORIENT_TILT),
    ("_rlnAnglePsi", md.RLN_ORIENT_PSI)
])
