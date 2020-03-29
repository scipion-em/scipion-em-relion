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

import os

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
import pwem
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtProcessParticles

import relion
import relion.convert as convert
from ..constants import STACK_NONE, STACK_MULT, STACK_ONE
from .protocol_base import ProtRelionBase


class ProtRelionExportParticles(ProtProcessParticles, ProtRelionBase):
    """ Export particles from Relion to be used outside Scipion. """

    _label = 'export particles'
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        form.addSection(label='Input')

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Select the input images from the project.')

        form.addParam('useAlignment', params.BooleanParam, default=True,
                      label='Write alignment information?',
                      help='If *Yes* the alignment information (2D or 3D) '
                           'will be written to the resulting .star file if '
                           'the particles contains such information.')

        form.addParam('stackType', params.EnumParam,
                      choices=["Don't write stacks",
                               "Write multiple stacks",
                               "Write a single stack"], default=STACK_MULT,
                      display=params.EnumParam.DISPLAY_LIST,
                      label="Binary stack files",
                      help="If *Don't write stacks* is chosen, only the star "
                           "files will be written out. Alternatively, you can "
                           "select to write images into a single stack file or"
                           " several stacks (one per micrograph). ")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        objId = self.inputParticles.get().getObjId()
        self._insertFunctionStep("exportParticlesStep", objId)

    # --------------------------- STEPS functions -----------------------------
    def exportParticlesStep(self, particlesId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        imgSet = self.inputParticles.get()
        self._stackType = self.stackType.get()
        self._ih = ImageHandler()
        self._stackDict = {}
        particlesPath = self._particlesPath()
        pwutils.cleanPath(particlesPath)
        pwutils.makePath(particlesPath)

        alignType = imgSet.getAlignment() if self.useAlignment else pwem.ALIGN_NONE
        outputDir = None
        outputStack = None
        postprocessImageRow = None

        if self._stackType == STACK_ONE:
            outputStack = self._particlesPath('particles.mrcs')
        elif self._stackType == STACK_MULT:
            outputDir = particlesPath

        # We still need to maintain some old code for Relion 3.0
        if relion.Plugin.IS_30():
            outputDir = self._getExtraPath()
            postprocessImageRow = self._postprocessImageRow

        # Create links to binary files and write the relion .star file
        convert.writeSetOfParticles(
            imgSet, self._particlesPath("particles.star"),
            outputDir=outputDir,
            outputStack=outputStack,
            alignType=alignType,
            postprocessImageRow=postprocessImageRow,
            fillMagnification=True,
            forceConvert=True)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []
        return validateMsgs
    
    def _summary(self):
        summary = []
        if pwutils.exists(self._particlesPath("particles.star")):
            summary.append('Particles were exported into: %s'
                           % self._particlesPath("particles.star"))
        else:
            summary.append('Output is not ready.')
        return summary
    
    # --------------------------- UTILS functions -----------------------------
    # TODO: Remove this function when support for 3.0 is dropped
    def _postprocessImageRow(self, img, row):
        """ Write the binary image to the final stack
        and update the row imageName. """

        if self._stackType > STACK_NONE:
            rlnImageName = row.getValue('rlnImageName')
            # backup the original name
            row.setValue('rlnOriginalParticleName', rlnImageName)

            if self._stackType == STACK_ONE:
                self._count = getattr(self, '_count', 1)
                index, stackName = (self._count, 'particles.mrcs')
                self._count += 1
            else:  # STACK_MULT
                baseName = pwutils.removeBaseExt(img.getFileName())
                if baseName not in self._stackDict:
                    self._stackDict[baseName] = 0
                index = self._stackDict[baseName] + 1
                stackName = baseName + '.mrcs'
                self._stackDict[baseName] = index

            stackFn = self._particlesPath(stackName)

            self._ih.convert(img, (index, stackFn))
            # Store relative path in the star file
            relStackFn = os.path.relpath(stackFn, self._getPath())
            row.setValue('rlnImageName',
                         convert.locationToRelion(index, relStackFn))

    def _particlesPath(self, *paths):
        return self._getPath('Particles', *paths)