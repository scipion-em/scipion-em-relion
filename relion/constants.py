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
RELION_CUDA_BIN = 'RELION_CUDA_BIN'
RELION_ENV_ACTIVATION = 'RELION_ENV_ACTIVATION'
DEFAULT_ACTIVATION_CMD = 'conda activate relion-%s'

SIDESPLITTER = 'SIDESPLITTER'
SIDESPLITTER_HOME = 'SIDESPLITTER_HOME'
RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE = 'RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE'
TORCH_HOME_VAR = 'TORCH_HOME'

# Supported versions:
V3_1 = '3.1'
V4_0 = '4.0'
V5_0 = '5.0'

MASK_FILL_ZERO = 0
MASK_FILL_NOISE = 1

ANGULAR_SAMPLING_LIST = ['30', '15', '7.5', '3.7', '1.8', '0.9', '0.5',
                         '0.2', '0.1', '0.06', '0.03', '0.01', '0.007', '0.004']

# Micrograph type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1

# Picking protocols
REF_AVERAGES = 0
REF_VOLUME = REF_BLOBS = 1

# Used in ctf-refinement
FIT_NO = 0
FIT_MICS = 1
FIT_PARTS = 2

# Axis code
AX_X = 0
AX_Y = 1
AX_Z = 2

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

FSC_TYPE_OVERALL = 0
FSC_TYPE_MODEL_MAP = 1
FSC_TYPE_WORK_FREE = 2

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

ALIGNMENT_DICT = OrderedDict([
    ("_rlnOriginX", md.RLN_ORIENT_ORIGIN_X),
    ("_rlnOriginY", md.RLN_ORIENT_ORIGIN_Y),
    ("_rlnOriginZ", md.RLN_ORIENT_ORIGIN_Z),
    ("_rlnAngleRot", md.RLN_ORIENT_ROT),
    ("_rlnAngleTilt", md.RLN_ORIENT_TILT),
    ("_rlnAnglePsi", md.RLN_ORIENT_PSI)
])


PARTICLE_EXTRA_LABELS = [
    'rlnCtfBfactor',
    'rlnCtfScalefactor',
    'rlnPhaseShift',
    'rlnGroupName',
    'rlnGroupNumber',
    'rlnRandomSubset'
]

# Labels from relion-4.0/src/metadata_label.h as of 20-Feb-2022

LABELS_DICT = {
    #"rlnAccumMotionEarly": float,
    #"rlnAccumMotionLate": float,
    #"rlnAccumMotionTotal": float,
    #"rlnAccuracyRotations": float,
    #"rlnAccuracyTranslations": float,
    #"rlnAccuracyTranslationsAngst": float,
    #"rlnAdaptiveOversampleFraction": float,
    #"rlnAdaptiveOversampleOrder": int,
    #"rlnAmplitudeContrast": float,
    #"rlnAmplitudeCorrelationMaskedMaps": float,
    #"rlnAmplitudeCorrelationUnmaskedMaps": float,
    #"rlnAnglePsi": float,
    #"rlnAnglePsiFlip": bool,
    #"rlnAnglePsiFlipRatio": float,
    #"rlnAnglePsiPrior": float,
    #"rlnAngleRot": float,
    #"rlnAngleRotFlipRatio": float,
    #"rlnAngleRotPrior": float,
    #"rlnAngleTilt": float,
    #"rlnAngleTiltPrior": float,
    #"rlnAngstromResolution": float,
    #"rlnAreaId": int,
    #"rlnAreaName": str,
    #"rlnAutoLocalSearchesHealpixOrder": int,
    #"rlnAutopickFigureOfMerit": float,
    #"rlnAvailableMemory": float,
    #"rlnAverageNrOfFrames": int,
    #"rlnAveragePmax": float,
    #"rlnAverageValue": float,
    #"rlnBeamTiltClass": int,
    #"rlnBeamTiltX": float,
    #"rlnBeamTiltY": float,
    #"rlnBestResolutionThusFar": float,
    #"rlnBfactorUsedForSharpening": float,
    #"rlnBodyKeepFixed": int,
    #"rlnBodyMaskName": str,
    #"rlnBodyReferenceName": str,
    #"rlnBodyRotateDirectionX": float,
    #"rlnBodyRotateDirectionY": float,
    #"rlnBodyRotateDirectionZ": float,
    #"rlnBodyRotateRelativeTo": int,
    #"rlnBodySigmaAngles": float,
    #"rlnBodySigmaOffset": float,
    #"rlnBodySigmaOffsetAngst": float,
    #"rlnBodySigmaPsi": float,
    #"rlnBodySigmaRot": float,
    #"rlnBodySigmaTilt": float,
    #"rlnBodyStarFile": str,
    #"rlnChangesOptimalClasses": float,
    #"rlnChangesOptimalOffsets": float,
    #"rlnChangesOptimalOrientations": float,
    #"rlnChromaticAberration": float,
    #"rlnCircleMaskedKurt": float,
    #"rlnCircleMaskedMean": float,
    #"rlnCircleMaskedSkew": float,
    #"rlnCircleMaskedStddev": float,
    #"rlnClassDistribution": float,
    #"rlnClassIndex": int,
    #"rlnClassNumber": int,
    #"rlnClassPriorOffsetX": float,
    #"rlnClassPriorOffsetY": float,
    #"rlnClassScore": float,
    #"rlnCoarseImageSize": int,
    #"rlnComment": str,
    #"rlnConvergenceCone": float,
    #"rlnCoordinateX": float,
    #"rlnCoordinateY": float,
    #"rlnCoordinateZ": float,
    #"rlnCorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps": float,
    #"rlnCorrelationFitGuinierPlot": float,
    #"rlnCtfAstigmatism": float,
    #"rlnCtfBfactor": float,
    #"rlnCtfDataAreCtfPremultiplied": bool,
    #"rlnCtfDataArePhaseFlipped": bool,
    #"rlnCtfFigureOfMerit": float,
    #"rlnCtfImage": str,
    #"rlnCtfMaxResolution": float,
    #"rlnCtfPowerSpectrum": str,
    #"rlnCtfScalefactor": float,
    #"rlnCtfValidationScore": float,
    #"rlnCtfValue": float,
    #"rlnCurrentImageSize": int,
    #"rlnCurrentIteration": int,
    #"rlnCurrentResolution": float,
    #"rlnDataDimensionality": int,
    #"rlnDataType": int,
    #"rlnDefocusAngle": float,
    #"rlnDefocusU": float,
    #"rlnDefocusV": float,
    #"rlnDetectorPixelSize": float,
    #"rlnDiff2RandomHalves": float,
    #"rlnDifferentialPhaseResidualMaskedMaps": float,
    #"rlnDifferentialPhaseResidualUnmaskedMaps": float,
    #"rlnDoAutoRefine": bool,
    #"rlnDoAutoSampling": bool,
    #"rlnDoCenterClasses": bool,
    #"rlnDoCorrectCtf": bool,
    #"rlnDoCorrectMagnification": bool,
    #"rlnDoCorrectNorm": bool,
    #"rlnDoCorrectScale": bool,
    #"rlnDoExternalReconstruct": bool,
    #"rlnDoFastSubsetOptimisation": bool,
    #"rlnDoGradientRefine": bool,
    #"rlnDoHelicalRefine": bool,
    #"rlnDoIgnoreCtfUntilFirstPeak": bool,
    #"rlnDoMapEstimation": bool,
    #"rlnDoOnlyFlipCtfPhases": bool,
    #"rlnDoRealignMovies": bool,
    #"rlnDoSkipAlign": bool,
    #"rlnDoSkipRotate": bool,
    #"rlnDoSolventFlattening": bool,
    #"rlnDoSolventFscCorrection": bool,
    #"rlnDoSplitRandomHalves": bool,
    #"rlnDoStochasticEM": bool,
    #"rlnDoStochasticGradientDescent": bool,
    #"rlnDoZeroMask": bool,
    #"rlnEERGrouping": int,
    #"rlnEERUpsampling": int,
    #"rlnEdgeSignal": float,
    #"rlnEnabled": bool,
    #"rlnEnergyLoss": float,
    #"rlnEstimatedResolution": float,
    #"rlnEvenZernike": str,
    #"rlnExperimentalDataStarFile": str,
    #"rlnExtReconsDataImag": str,
    #"rlnExtReconsDataReal": str,
    #"rlnExtReconsResult": str,
    #"rlnExtReconsResultStarfile": str,
    #"rlnExtReconsWeight": str,
    #"rlnFftKurt": float,
    #"rlnFftMean": float,
    #"rlnFftSkew": float,
    #"rlnFftStddev": float,
    #"rlnFinalResolution": float,
    #"rlnFittedInterceptGuinierPlot": float,
    #"rlnFittedSlopeGuinierPlot": float,
    #"rlnFixSigmaNoiseEstimates": bool,
    #"rlnFixSigmaOffsetEstimates": bool,
    #"rlnFixTauEstimates": bool,
    #"rlnFourierCompleteness": float,
    #"rlnFourierMask": str,
    #"rlnFourierShellCorrelation": float,
    #"rlnFourierShellCorrelationCorrected": float,
    #"rlnFourierShellCorrelationMaskedMaps": float,
    #"rlnFourierShellCorrelationParticleMaskFraction": float,
    #"rlnFourierShellCorrelationParticleMolWeight": float,
    #"rlnFourierShellCorrelationUnmaskedMaps": float,
    #"rlnFourierSpaceInterpolator": int,
    #"rlnGoldStandardFsc": float,
    #"rlnGradCurrentStepsize": float,
    #"rlnGradEmIters": int,
    #"rlnGradHasConverged": bool,
    #"rlnGradMoment1": str,
    #"rlnGradMoment2": str,
    #"rlnGradSubsetOrder": int,
    #"rlnGradSuspendFinerSamplingIter": int,
    #"rlnGradSuspendLocalSamplingIter": int,
    #"rlnGranulo": str,
    "rlnGroupName": str,
    #"rlnGroupNrParticles": int,
    "rlnGroupNumber": int,
    #"rlnGroupScaleCorrection": float,
    #"rlnHasConverged": bool,
    #"rlnHasHighFscAtResolLimit": bool,
    #"rlnHasLargeSizeIncreaseIterationsAgo": int,
    #"rlnHealpixOrder": int,
    #"rlnHealpixOrderOriginal": int,
    #"rlnHelicalCentralProportion": float,
    #"rlnHelicalKeepTiltPriorFixed": bool,
    #"rlnHelicalMaskTubeInnerDiameter": float,
    #"rlnHelicalMaskTubeOuterDiameter": float,
    #"rlnHelicalOffsetStep": float,
    #"rlnHelicalRise": float,
    #"rlnHelicalRiseInitial": float,
    #"rlnHelicalRiseInitialStep": float,
    #"rlnHelicalRiseMax": float,
    #"rlnHelicalRiseMin": float,
    #"rlnHelicalSigmaDistance": float,
    #"rlnHelicalSymmetryLocalRefinement": bool,
    #"rlnHelicalTrackLength": float,
    #"rlnHelicalTrackLengthAngst": float,
    #"rlnHelicalTubeID": int,
    #"rlnHelicalTubePitch": float,
    #"rlnHelicalTwist": float,
    #"rlnHelicalTwistInitial": float,
    #"rlnHelicalTwistInitialStep": float,
    #"rlnHelicalTwistMax": float,
    #"rlnHelicalTwistMin": float,
    #"rlnHighresLimitExpectation": float,
    #"rlnHighresLimitSGD": float,
    #"rlnIgnoreHelicalSymmetry": bool,
    #"rlnImageDimensionality": int,
    #"rlnImageId": int,
    #"rlnImageName": str,
    #"rlnImageOriginalName": str,
    #"rlnImagePixelSize": float,
    #"rlnImageSize": int,
    #"rlnImageSizeX": int,
    #"rlnImageSizeY": int,
    #"rlnImageSizeZ": int,
    #"rlnImageWeight": float,
    #"rlnIncrementImageSize": int,
    #"rlnInnerCircleKurt": float,
    #"rlnInnerCircleMean": float,
    #"rlnInnerCircleSkew": float,
    #"rlnInnerCircleStddev": float,
    #"rlnIs3DSampling": bool,
    #"rlnIs3DTranslationalSampling": bool,
    #"rlnIsFlip": bool,
    #"rlnIsHelix": bool,
    #"rlnIsSelected": int,
    #"rlnJobIsContinue": bool,
    #"rlnJobIsTomo": bool,
    #"rlnJobOptionDefaultValue": str,
    #"rlnJobOptionDirectoryDefault": str,
    #"rlnJobOptionFilePattern": str,
    #"rlnJobOptionGUILabel": str,
    #"rlnJobOptionHelpText": str,
    #"rlnJobOptionMenuOptions": str,
    #"rlnJobOptionSliderMax": float,
    #"rlnJobOptionSliderMin": float,
    #"rlnJobOptionSliderStep": float,
    #"rlnJobOptionValue": str,
    #"rlnJobOptionVariable": str,
    #"rlnJobScore": float,
    #"rlnJobType": int,
    #"rlnJobTypeLabel": str,
    #"rlnJoboptionType": int,
    #"rlnJoinHalvesUntilThisResolution": float,
    #"rlnKullbackLeiblerDivergence": float,
    #"rlnKullbackLeibnerDivergence": float,
    #"rlnKurtosisExcessValue": float,
    #"rlnLBP": str,
    #"rlnLensStability": float,
    #"rlnLocalSymmetryFile": str,
    #"rlnLogAmplitudesIntercept": float,
    #"rlnLogAmplitudesMTFCorrected": float,
    #"rlnLogAmplitudesOriginal": float,
    #"rlnLogAmplitudesSharpened": float,
    #"rlnLogAmplitudesWeighted": float,
    #"rlnLogLikeliContribution": float,
    #"rlnLogLikelihood": float,
    #"rlnLongitudinalDisplacement": float,
    #"rlnLowpassFilteredImageMax": float,
    #"rlnLowpassFilteredImageMean": float,
    #"rlnLowpassFilteredImageMin": float,
    #"rlnLowpassFilteredImageStddev": float,
    #"rlnLowresLimitExpectation": float,
    #"rlnMagMat00": float,
    #"rlnMagMat01": float,
    #"rlnMagMat10": float,
    #"rlnMagMat11": float,
    #"rlnMagnification": float,
    #"rlnMagnificationCorrection": float,
    #"rlnMagnificationSearchRange": float,
    #"rlnMagnificationSearchStep": float,
    #"rlnMaskName": str,
    #"rlnMatrix_1_1": float,
    #"rlnMatrix_1_2": float,
    #"rlnMatrix_1_3": float,
    #"rlnMatrix_2_1": float,
    #"rlnMatrix_2_2": float,
    #"rlnMatrix_2_3": float,
    #"rlnMatrix_3_1": float,
    #"rlnMatrix_3_2": float,
    #"rlnMatrix_3_3": float,
    #"rlnMaxNumberOfPooledParticles": int,
    #"rlnMaxValueProbDistribution": float,
    #"rlnMaximumCoarseImageSize": int,
    #"rlnMaximumValue": float,
    #"rlnMicrographBinning": float,
    #"rlnMicrographCoordinates": str,
    #"rlnMicrographDefectFile": str,
    #"rlnMicrographDoseRate": float,
    #"rlnMicrographEndFrame": int,
    #"rlnMicrographFrameNumber": int,
    #"rlnMicrographGainName": str,
    #"rlnMicrographId": int,
    #"rlnMicrographMetadata": str,
    "rlnMicrographMovieName": str,
    "rlnMicrographName": str,
    "rlnMicrographNameNoDW": str,
    #"rlnMicrographOriginalPixelSize": float,
    #"rlnMicrographPixelSize": float,
    #"rlnMicrographPreExposure": float,
    #"rlnMicrographShiftX": float,
    #"rlnMicrographShiftY": float,
    #"rlnMicrographStartFrame": int,
    #"rlnMicrographTiltAngle": float,
    #"rlnMicrographTiltAxisDirection": float,
    #"rlnMicrographTiltAxisOutOfPlane": float,
    #"rlnMinRadiusNnInterpolation": int,
    #"rlnMinimumValue": float,
    #"rlnModelStarFile": str,
    #"rlnModelStarFile2": str,
    #"rlnMolecularWeight": float,
    #"rlnMotionModelCoeff": float,
    #"rlnMotionModelCoeffsIdx": int,
    #"rlnMotionModelVersion": int,
    #"rlnMovieFrameNumber": int,
    #"rlnMovieFramesRunningAverage": int,
    #"rlnMtfFileName": str,
    #"rlnMtfValue": float,
    #"rlnNormCorrection": float,
    #"rlnNormCorrectionAverage": float,
    #"rlnNormalizedFeatureVector": str,
    #"rlnNrBodies": int,
    #"rlnNrClasses": int,
    #"rlnNrGroups": int,
    #"rlnNrHelicalAsymUnits": int,
    #"rlnNrHelicalNStart": int,
    #"rlnNrOfFrames": int,
    #"rlnNrOfSignificantSamples": int,
    #"rlnNrOpticsGroups": int,
    #"rlnNumberOfIterWithoutChangingAssignments": int,
    #"rlnNumberOfIterWithoutResolutionGain": int,
    #"rlnNumberOfIterations": int,
    #"rlnOddZernike": str,
    #"rlnOffsetRange": float,
    #"rlnOffsetRangeOriginal": float,
    #"rlnOffsetStep": float,
    #"rlnOffsetStepOriginal": float,
    "rlnOpticsGroup": int,
    "rlnOpticsGroupName": str,
    #"rlnOpticsGroupNrParticles": int,
    #"rlnOpticsGroupNumber": int,
    #"rlnOpticsStarFile": str,
    #"rlnOrientSamplingStarFile": str,
    #"rlnOrientationDistribution": float,
    #"rlnOrientationalPriorMode": int,
    #"rlnOrientationsID": int,
    #"rlnOriginX": float,
    #"rlnOriginXAngst": float,
    #"rlnOriginXPrior": float,
    #"rlnOriginXPriorAngst": float,
    #"rlnOriginY": float,
    #"rlnOriginYAngst": float,
    #"rlnOriginYPrior": float,
    #"rlnOriginYPriorAngst": float,
    #"rlnOriginZ": float,
    #"rlnOriginZAngst": float,
    #"rlnOriginZPrior": float,
    #"rlnOriginZPriorAngst": float,
    #"rlnOriginalImageSize": int,
    #"rlnOriginalParticleName": str,
    #"rlnOutputRootName": str,
    #"rlnOverallAccuracyRotations": float,
    #"rlnOverallAccuracyTranslations": float,
    #"rlnOverallAccuracyTranslationsAngst": float,
    #"rlnOverallFourierCompleteness": float,
    #"rlnPaddingFactor": float,
    #"rlnParticleBoxFractionMolecularWeight": float,
    #"rlnParticleBoxFractionSolventMask": float,
    #"rlnParticleDiameter": float,
    #"rlnParticleFigureOfMerit": float,
    #"rlnParticleId": int,
    #"rlnParticleName": str,
    #"rlnParticleNr": float,
    #"rlnParticleNumber": int,
    #"rlnParticleSelectZScore": float,
    #"rlnPerFrameCumulativeWeight": float,
    #"rlnPerFrameRelativeWeight": float,
    #"rlnPhaseShift": float,
    #"rlnPipeLineEdgeFromNode": str,
    #"rlnPipeLineEdgeProcess": str,
    #"rlnPipeLineEdgeToNode": str,
    #"rlnPipeLineJobCounter": int,
    #"rlnPipeLineNodeName": str,
    #"rlnPipeLineNodeType": int,
    #"rlnPipeLineNodeTypeLabel": str,
    #"rlnPipeLineProcessAlias": str,
    #"rlnPipeLineProcessName": str,
    #"rlnPipeLineProcessStatus": int,
    #"rlnPipeLineProcessStatusLabel": str,
    #"rlnPipeLineProcessType": int,
    #"rlnPipeLineProcessTypeLabel": str,
    #"rlnPixelSize": float,
    #"rlnPostprocessedMap": str,
    #"rlnPostprocessedMapMasked": str,
    #"rlnPredictedClassScore": float,
    #"rlnProteinArea": int,
    #"rlnProteinCAR": float,
    #"rlnProteinEntropy": float,
    #"rlnProteinHaralick": str,
    #"rlnProteinKurt": float,
    #"rlnProteinLBP": str,
    #"rlnProteinMean": float,
    #"rlnProteinSkew": float,
    #"rlnProteinStddev": float,
    #"rlnProteinSum": float,
    #"rlnPsiStep": float,
    #"rlnPsiStepOriginal": float,
    #"rlnRadiusMaskExpImages": int,
    #"rlnRadiusMaskMap": int,
    #"rlnRandomSeed": int,
    #"rlnRandomSubset": int,
    #"rlnRandomiseFrom": float,
    #"rlnReconstructImageName": str,
    #"rlnReferenceDimensionality": int,
    #"rlnReferenceImage": str,
    #"rlnReferenceSigma2": float,
    #"rlnReferenceSpectralPower": float,
    #"rlnReferenceTau2": float,
    #"rlnRefsAreCtfCorrected": bool,
    #"rlnRelativeResolution": float,
    #"rlnRelativeSignalIntensity": float,
    #"rlnResolution": float,
    #"rlnResolutionInversePixel": float,
    #"rlnResolutionSquared": float,
    #"rlnRingKurt": float,
    #"rlnRingMean": float,
    #"rlnRingSkew": float,
    #"rlnRingStddev": float,
    #"rlnSamplingPerturbFactor": float,
    #"rlnSamplingPerturbInstance": float,
    #"rlnSamplingRate": float,
    #"rlnSamplingRateX": float,
    #"rlnSamplingRateY": float,
    #"rlnSamplingRateZ": float,
    #"rlnScatteredSignal": float,
    #"rlnSchemeBooleanVariableName": str,
    #"rlnSchemeBooleanVariableResetValue": bool,
    #"rlnSchemeBooleanVariableValue": bool,
    #"rlnSchemeCurrentNodeName": str,
    #"rlnSchemeEdgeBooleanVariable": str,
    #"rlnSchemeEdgeInputNodeName": str,
    #"rlnSchemeEdgeIsFork": bool,
    #"rlnSchemeEdgeNumber": int,
    #"rlnSchemeEdgeOutputNodeName": str,
    #"rlnSchemeEdgeOutputNodeNameIfTrue": str,
    #"rlnSchemeEmailAddress": str,
    #"rlnSchemeFloatVariableName": str,
    #"rlnSchemeFloatVariableResetValue": float,
    #"rlnSchemeFloatVariableValue": float,
    #"rlnSchemeJobHasStarted": bool,
    #"rlnSchemeJobMode": str,
    #"rlnSchemeJobName": str,
    #"rlnSchemeJobNameOriginal": str,
    #"rlnSchemeName": str,
    #"rlnSchemeOperatorInput1": str,
    #"rlnSchemeOperatorInput2": str,
    #"rlnSchemeOperatorName": str,
    #"rlnSchemeOperatorOutput": str,
    #"rlnSchemeOperatorType": str,
    #"rlnSchemeOriginalStartNodeName": str,
    #"rlnSchemeStringVariableName": str,
    #"rlnSchemeStringVariableResetValue": str,
    #"rlnSchemeStringVariableValue": str,
    #"rlnSelected": int,
    #"rlnSgdClassInactivityThreshold": float,
    #"rlnSgdFinalIterationsFraction": float,
    #"rlnSgdFinalResolution": float,
    #"rlnSgdFinalSubsetSize": int,
    #"rlnSgdInitialIterationsFraction": float,
    #"rlnSgdInitialResolution": float,
    #"rlnSgdInitialSubsetSize": int,
    #"rlnSgdMaxSubsets": int,
    #"rlnSgdMinimumResolution": float,
    #"rlnSgdMuFactor": float,
    #"rlnSgdSigma2FudgeHalflife": int,
    #"rlnSgdSigma2FudgeInitial": float,
    #"rlnSgdSkipAnneal": bool,
    #"rlnSgdStepsize": float,
    "rlnSgdStepsizeScheme": str,  # used by two labels, str and int
    #"rlnSgdSubsetSize": int,
    #"rlnSgdWriteEverySubset": int,
    #"rlnSigma2Noise": float,
    #"rlnSigmaOffsets": float,
    #"rlnSigmaOffsetsAngst": float,
    #"rlnSigmaPriorPsiAngle": float,
    #"rlnSigmaPriorRotAngle": float,
    #"rlnSigmaPriorTiltAngle": float,
    #"rlnSignalToNoiseRatio": float,
    #"rlnSkewnessValue": float,
    #"rlnSmallestChangesClasses": int,
    #"rlnSmallestChangesOffsets": float,
    #"rlnSmallestChangesOrientations": float,
    #"rlnSolventArea": int,
    #"rlnSolventEntropy": float,
    #"rlnSolventHaralick": str,
    #"rlnSolventKurt": float,
    #"rlnSolventLBP": str,
    #"rlnSolventMask2Name": str,
    #"rlnSolventMaskName": str,
    #"rlnSolventMean": float,
    #"rlnSolventSkew": float,
    #"rlnSolventStddev": float,
    #"rlnSolventSum": float,
    #"rlnSortedIndex": int,
    #"rlnSpectralIndex": int,
    #"rlnSpectralOrientabilityContribution": float,
    #"rlnSphericalAberration": float,
    #"rlnSsnrMap": float,
    #"rlnStandardDeviationValue": float,
    #"rlnStarFileMovieParticles": str,
    #"rlnSubImageStack": str,
    #"rlnSubImageStarFile": str,
    #"rlnSymmetryGroup": str,
    #"rlnTau2FudgeArg": float,
    #"rlnTau2FudgeFactor": float,
    #"rlnTau2FudgeScheme": str,
    #"rlnTauSpectrumName": str,
    #"rlnTiltAngleLimit": float,
    #"rlnTomoDefocusSlope": float,
    #"rlnTomoDeformationCoefficients": str,
    #"rlnTomoDeformationGridSizeX": int,
    #"rlnTomoDeformationGridSizeY": int,
    #"rlnTomoDeformationType": str,
    #"rlnTomoFiducialsStarFile": str,
    #"rlnTomoFrameCount": int,
    #"rlnTomoHand": float,
    #"rlnTomoIceNormalX": float,
    #"rlnTomoIceNormalY": float,
    #"rlnTomoIceNormalZ": float,
    #"rlnTomoImportCtfFindFile": str,
    #"rlnTomoImportCtfPlotterFile": str,
    #"rlnTomoImportCulledFile": str,
    #"rlnTomoImportFractionalDose": float,
    #"rlnTomoImportImodDir": str,
    #"rlnTomoImportOffsetX": float,
    #"rlnTomoImportOffsetY": float,
    #"rlnTomoImportOffsetZ": float,
    #"rlnTomoImportOrderList": str,
    #"rlnTomoImportParticleFile": str,
    #"rlnTomoManifoldIndex": int,
    #"rlnTomoManifoldParams": str,
    #"rlnTomoManifoldType": str,
    #"rlnTomoManifoldsFile": str,
    "rlnTomoName": str,
    #"rlnTomoParticleId": int,
    #"rlnTomoParticleName": str,
    #"rlnTomoParticlesFile": str,
    #"rlnTomoProjW": str,
    #"rlnTomoProjX": str,
    #"rlnTomoProjY": str,
    #"rlnTomoProjZ": str,
    #"rlnTomoReferenceFscFile": str,
    #"rlnTomoReferenceMap1File": str,
    #"rlnTomoReferenceMap2File": str,
    #"rlnTomoReferenceMaskFile": str,
    #"rlnTomoRelativeIceThickness": float,
    #"rlnTomoRelativeLuminance": float,
    #"rlnTomoSizeX": int,
    #"rlnTomoSizeY": int,
    #"rlnTomoSizeZ": int,
    #"rlnTomoSubtomogramBinning": float,
    #"rlnTomoSubtomogramPsi": float,
    #"rlnTomoSubtomogramRot": float,
    #"rlnTomoSubtomogramTilt": float,
    #"rlnTomoTempPredSquared": float,
    #"rlnTomoTempPredTimesObs": float,
    #"rlnTomoTiltMovieFile": str,
    #"rlnTomoTiltMovieIndex": int,
    "rlnTomoTiltSeriesName": str,
    #"rlnTomoTiltSeriesPixelSize": float,
    #"rlnTomoTomogramsFile": str,
    #"rlnTomoTrajectoriesFile": str,
    #"rlnTotalEntropy": float,
    #"rlnTransversalDisplacement": float,
    #"rlnUnfilteredMapHalf1": str,
    #"rlnUnfilteredMapHalf2": str,
    #"rlnUnknownLabel": str,
    #"rlnUseTooCoarseSampling": bool,
    #"rlnVoltage": float,
    #"rlnWeightedResolution": float,
    #"rlnWidthMaskEdge": int,
    #"rlnZernikeMoments": str,
}
