# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
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

import relion.convert as convert
from .protocol_base import ProtRelionBase


class ProtRelionExportCoordinates(ProtRelionBase):
    """ Export coordinates from Relion to be used outside Scipion. """

    _label = 'export coordinates'
    _devStatus = PROD

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label="Input coordinates",
                      help='Select the SetOfCoordinates ')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        objId = self.inputCoordinates.get().getObjId()
        self._insertFunctionStep("exportCoordsStep", objId)

    # --------------------------- STEPS functions -----------------------------
    def exportCoordsStep(self, coordsId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        """
        pwutils.cleanPath(self._getExportPath())
        pwutils.makePath(self._getExportPath())

        convert.writeSetOfCoordinates(self._getExportPath(),
                                      self.getCoords(),
                                      self._getMicPos)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        return validateMsgs

    def _summary(self):
        summary = []

        summary.append("Output is written to: \n%s\n" %
                       os.path.abspath(self._getExportPath()))

        return summary
    
    # --------------------------- UTILS functions -----------------------------
    def _getMicPos(self, mic):
        fileName = pwutils.removeBaseExt(mic.getFileName()) + "_coords.star"
        return fileName

    def getCoords(self):
        return self.inputCoordinates.get()

    def _getExportPath(self, *paths):
        return os.path.join(self._getPath('Export'), *paths)
