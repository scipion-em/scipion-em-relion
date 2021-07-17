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
from pyworkflow.constants import PROD
from pwem.constants import ALIGN_NONE
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtProcessParticles

import relion.convert as convert
from ..constants import STACK_MULT, STACK_ONE
from .protocol_base import ProtRelionBase


class ProtRelionExportParticles(ProtProcessParticles, ProtRelionBase):
    """ Export particles from Relion to be used outside Scipion. """

    _label = 'export particles'
    _devStatus = PROD
    PTCLS_STAR_FILE = 'particles_%06d.star'
    
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
        pwutils.cleanPath(self._getExportPath())
        pwutils.makePath(self._getExportPath())
        imgSet = self.inputParticles.get()
        self._stackType = self.stackType.get()
        self._ih = ImageHandler()
        self._stackDict = {}

        alignType = imgSet.getAlignment() if self.useAlignment else ALIGN_NONE
        outputDir = None
        outputStack = None
        postprocessImageRow = None

        if self._stackType == STACK_ONE:
            outputStack = self._getExportPath('Particles/particles.mrcs')
            pwutils.makePath(self._getExportPath("Particles"))

        elif self._stackType == STACK_MULT:
            postprocessImageRow = self._postprocessImageRow
            outputDir = self._getExportPath("Particles")

        # Create links to binary files and write the relion .star file
        convert.writeSetOfParticles(
            imgSet, self._getStarFile(),
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

        if os.path.exists(self._getStarFile()):
            summary.append("Output is written to: \n%s\n" %
                           os.path.abspath(self._getExportPath()))
            summary.append("Pixel size: *%0.3f*" % self._getPixelSize())
        else:
            summary.append("No output generated yet.")

        return summary
    
    # --------------------------- UTILS functions -----------------------------
    def _postprocessImageRow(self, img, row):
        """ Stack fn should be relative to Export.
        Only relevant when saving multiple stacks. """
        convert.relativeFromFileName(row, self._getExportPath())

    def _getExportPath(self, *paths):
        return os.path.join(self._getPath('Export'), *paths)

    def _getStarFile(self):
        return self._getExportPath(self.PTCLS_STAR_FILE % self.getObjId())

    def _getPixelSize(self):
        return self.inputParticles.get().getSamplingRate()
