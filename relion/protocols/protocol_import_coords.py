# **************************************************************************
# *
# * Authors: Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from enum import Enum

from pwem.protocols import EMProtocol
from pwem.objects import SetOfCoordinates
from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils

from relion import convert


class outputs(Enum):
    outputCoordinates = SetOfCoordinates


class ProtRelionImportCoords(EMProtocol):
    """
    Import coordinates from a particles star file.
    """
    _label = 'import coordinates'
    _devStatus = BETA
    _possibleOutputs = outputs

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Import')
        form.addParam('filePath', params.FileParam,
                      allowsNull=False,
                      label="Input STAR file path with coordinates")
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label='Input micrographs',
                      help='Select the micrographs for which you want to '
                           'import coordinates.')
        form.addParam('boxSize', params.IntParam, label='Box size')
        form.addParam('scale', params.FloatParam,
                      label='Scale', default=1,
                      help='Factor to scale coordinates')
        form.addParam('invertX', params.BooleanParam, default=False,
                      label='Invert X')
        form.addParam('invertY', params.BooleanParam, default=False,
                      label='Invert Y',
                      help='Invert Y for EMAN coordinates taken on dm3 or'
                           ' tif micrographs')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep', self.filePath.get())

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self, coordStar):
        inputMics = self.getMicrographList()
        micList = [pwutils.removeExt(mic.getMicName()) for mic in inputMics]
        coordsSet = self._createSetOfCoordinates(inputMics)
        coordsSet.setBoxSize(self.boxSize.get())

        reader = convert.createReader()
        reader.readSetOfCoordinates(
            coordStar, coordsSet, micList,
            postprocessCoordRow=self._postprocessCoordRow)

        self._defineOutputs(**{outputs.outputCoordinates.name: coordsSet})
        self._defineSourceRelation(self.inputMicrographs, coordsSet)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        if not hasattr(self, 'outputCoordinates'):
            msg = 'Output coordinates not ready yet'
        else:
            msg = "%d coordinates from " % self.outputCoordinates.getSize()
            msg += "micrographs %s " % self.getObjectTag('inputMicrographs')
            if self.scale.get() != 1.:
                msg += " Scale factor %0.2f was applied." % self.scale
            if self.invertX.get():
                msg += " X coordinate was inverted."
            if self.invertY.get():
                msg += " Y coordinate was inverted."

            msg += "Output coordinates: %s." % self.getObjectTag('outputCoordinates')

        return [msg]

    # -------------------------- UTILS functions ------------------------------
    def _postprocessCoordRow(self, coord, row):
        scaleFactor = self.scale.get()
        x = coord.getX()
        y = coord.getY()
        if scaleFactor != 1.:
            x = coord.getX() * scaleFactor
            y = coord.getY() * scaleFactor
        if self.invertX:
            width = self.getMicDim()[0]
            x = width - x
        if self.invertY:
            height = self.getMicDim()[1]
            y = height - y
        coord.setX(x)
        coord.setY(y)

    def getMicrographList(self):
        return self.inputMicrographs.get()

    def getMicDim(self):
        return self.getMicrographList().getDimensions()
