# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow.utils as pwutils
from pyworkflow.protocol.params import (PointerParam, BooleanParam,
                                        FloatParam, IntParam, Positive)
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.constants import PROD
from pwem.constants import NO_INDEX
from pwem.objects import SetOfAverages
from pwem.protocols import ProtProcessParticles

import relion.convert as convert
from .protocol_base import ProtRelionBase


class ProtRelionPreprocessParticles(ProtProcessParticles, ProtRelionBase):
    """ This protocol wraps relion_preprocess program.

    It is used to perform normalisation, filtering or scaling of
    the particles.
    """
    _label = 'preprocess particles'
    _devStatus = PROD
    
    def __init__(self, **args):
        ProtProcessParticles.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Select the input images from the project.')
        
        form.addParam('doNormalize', BooleanParam, default=True,
                      label='Normalize', important=True,
                      help='If set to True, particles will be normalized in the'
                           'way RELION prefers it. It is recommended to '
                           '*always normalize your particles*, and use a '
                           'reasonable radius for the circle around your '
                           'particles outside of which the standard deviation '
                           'and average values for the noise are calculated.\n'
                           '*Note*: if the particles are re-scaled, the radius '
                           'for normalize will be taken over the new '
                           'dimensions.')
        form.addParam('backRadius', IntParam, default=-1,
                      condition='doNormalize',
                      label='Background radius (px)',
                      help='Pixels outside this circle are assumed to be '
                           'noise and their stddev is set to 1. Radius for '
                           'background circle definition (in pixel).')
        
        form.addParam('doRemoveDust', BooleanParam, default=False,
                      label='Remove dust from particles',
                      help='If there are white or black artifacts on the '
                           'micrographs (e.g. caused by dust or hot/dead '
                           'pixels), these may be removed by using a positive '
                           'value for the dust removal options. All '
                           'black/white pixels with values above the given '
                           'parameter times the standard deviation of the '
                           'noise are replaced by random values from a'
                           'Gaussian distribution. For cryo-EM data, values'
                           'around 3.5-5 are often useful. Make sure you do '
                           'not erase part of the true signal.')
        line = form.addLine('Dust sigmas', condition='doRemoveDust',
                            help='Sigma-values above which white/black dust '
                                 'will be removed (negative value means no '
                                 'dust removal)')
        line.addParam('whiteDust', FloatParam, default=-1., label='White')
        line.addParam('blackDust', FloatParam, default=-1., label='Black')
        
        form.addParam('doInvert', BooleanParam, default=False,
                      label='Invert contrast',
                      help='Invert the contrast if your particles are black '
                           'over a white background.')
        
        form.addSection('Scale and window')
        form.addParam('doScale', BooleanParam, default=False,
                      label='Scale particles?',
                      help='Re-scale the particles to this size (in pixels).')
        form.addParam('scaleSize', IntParam, default=0, validators=[Positive],
                      condition='doScale',
                      label='Scale size (px)',
                      help='New particle size in pixels.')
        form.addParam('doWindow', BooleanParam, default=False,
                      label='Window particles?',
                      help='Re-window the particles to this size (in pixels).')
        form.addParam('windowSize', IntParam, default=0, validators=[Positive],
                      condition='doWindow',
                      label='Window size (px)',
                      help='New particles windows size (in pixels).')

        form.addParallelSection(threads=4, mpi=1)
    
    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        inputParts = self.inputParticles.get()
        objId = self.inputParticles.get().getObjId()
        stackFiles = list(inputParts.getFiles())
        allIds = []
        args = self._getArgs()

        for stack in sorted(stackFiles):
            allIds.append(self._insertFunctionStep('processStep', objId,
                                                   stack, args,
                                                   prerequisites=[]))

        self._insertFunctionStep('createOutputStep', prerequisites=allIds)

    # --------------------------- STEPS functions -----------------------------
    def _getArgs(self):
        args = ''

        if self.doNormalize:
            radius = self.backRadius.get()
            if radius <= 0:
                radius = self._getOutputRadius()
            args += ' --norm --bg_radius %d' % radius

        if self.doRemoveDust:
            wDust = self.whiteDust.get()
            if wDust > 0:
                args += ' --white_dust %f' % wDust
            bDust = self.blackDust.get()
            if bDust > 0:
                args += ' --black_dust %f' % bDust

        if self.doInvert:
            args += ' --invert_contrast'

        if self.doScale:
            args += ' --scale %d' % self.scaleSize

        if self.doWindow:
            args += ' --window %d' % self.windowSize

        return args

    def processStep(self, objId, stack, args):
        # Enter here to generate the star file or to preprocess the images
        stackOut = self._getOutStack(stack)
        params = '--operate_on %s %s --operate_out %s' % (stack, args, stackOut)
        self.runJob(self._getProgram('relion_preprocess'), params)
        
    def createOutputStep(self):
        inputSet = self.inputParticles.get()
        
        if isinstance(inputSet, SetOfAverages):
            imgSet = self._createSetOfAverages()
        else:
            imgSet = self._createSetOfParticles()
        
        imgSet.copyInfo(inputSet)

        if self.doScale:
            oldSampling = inputSet.getSamplingRate()
            scaleFactor = self._getScaleFactor(inputSet)
            newSampling = oldSampling * scaleFactor
            imgSet.setSamplingRate(newSampling)

        imgSet.copyItems(inputSet, updateItemCallback=self._setFileName)
        self._defineOutputs(outputParticles=imgSet)
        self._defineTransformRelation(inputSet, imgSet)
    
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        if self.doScale and self.scaleSize.get() % 2 != 0:
            validateMsgs.append("Only re-scaling to even-sized images is "
                                "allowed in RELION.")
        
        if self.doWindow and self.windowSize.get() % 2 != 0:
            validateMsgs.append("Only re-windowing to even-sized images is "
                                "allowed in RELION.")
        
        if self.doNormalize:
            outputRadius = self._getOutputRadius()
            if self.backRadius > outputRadius:
                validateMsgs.append('Set a normalization background radius '
                                    'less than the particles output radius '
                                    '(%d pixels).' % outputRadius)
        return validateMsgs
    
    def _summary(self):
        summary = list()
        summary.append('Operations applied:')
        if self.doScale:
            summary.append(
                "- Particles scaled to *%d* pixels" % self.scaleSize.get())
        if self.doWindow:
            summary.append(
                "- Particles windowed to *%d* pixels" % self.windowSize.get())
        if self.doNormalize:
            summary.append("- Normalization, background radius *%d* pixels"
                           % self.backRadius.get())
        if self.doRemoveDust:
            if self.whiteDust > 0:
                summary.append(
                    '- Removed white dust (sigma=%0.3f)' % self.whiteDust.get())
            if self.blackDust > 0:
                summary.append(
                    '- Removed black dust (sigma=%0.3f)' % self.blackDust.get())
        if self.doInvert:
            summary.append('- Inverted contrast')
        
        return summary
    
    def _methods(self):
        summary = 'The input particles were preprocessed using Relion.'
        if self.doScale:
            summary += " Scaled to *%d* pixels." % self.scaleSize.get()
        if self.doWindow:
            summary += " Windowed to *%d* pixels." % self.windowSize.get()
        if self.doNormalize:
            summary += (" The particles were normalized using a background "
                        "radius of *%d* pixels." % self.backRadius.get())
        if self.doRemoveDust:
            if self.whiteDust > 0:
                summary += (' White dust was removed (pixels with a sigma > '
                            '%0.3f).' % self.whiteDust.get())
            if self.blackDust > 0:
                summary += (' Black dust was removed (pixels with a sigma > '
                            '%0.3f).' % self.blackDust.get())
        if self.doInvert:
            summary += ' The original constrast was inverted.'
        
        return [summary]
    # --------------------------- UTILS functions -----------------------------
    
    def _getOutputRadius(self):
        """ Get the radius of the output particles"""
        if self.doScale:
            radius = self.scaleSize.get() / 2
        else:
            xdim = self.inputParticles.get().getDimensions()[0]
            radius = xdim / 2
        return radius
    
    def _postprocessImageRow(self, img, imgRow):
        """ Since relion_preprocess will runs in its working directory
        we need to modify the default image path (from project dir)
        and make them relative to run working dir.
        """
        convert.relativeFromFileName(imgRow, self._getPath())
    
    def _getScaleFactor(self, inputSet):
        xdim = self.inputParticles.get().getDim()[0]
        scaleFactor = xdim / float(
            self.scaleSize.get() if self.doScale else xdim)
        return scaleFactor

    def _getOutStack(self, stack):
        """ Return the output stack filename based on the input. """
        return self._getExtraPath(pwutils.replaceBaseExt(stack, 'mrcs'))

    def _setFileName(self, item, row=None):
        index, fn = item.getLocation()
        index = 1 if index == NO_INDEX else index
        item.setLocation(index, self._getOutStack(fn))
        
        invFactor = 1 / self._getScaleFactor(item)
        
        if invFactor != 1.0:
            if item.hasCoordinate():
                item.scaleCoordinate(invFactor)
            if item.hasTransform():
                item.getTransform().scaleShifts(invFactor)
