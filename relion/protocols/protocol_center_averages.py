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

from enum import Enum

from pwem.protocols import ProtProcessParticles
from pwem.objects import SetOfAverages
from pyworkflow.protocol.params import PointerParam
from pyworkflow.constants import PROD

from .protocol_base import ProtRelionBase


class outputs(Enum):
    outputAverages = SetOfAverages


class ProtRelionCenterAverages(ProtProcessParticles, ProtRelionBase):
    """
    Align class averages by their center of mass using *relion_image_handler*.
     (With *--shift_com* option)
    """
    _label = 'center averages'
    _devStatus = PROD
    _possibleOutputs = outputs
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputAverages', PointerParam,
                      pointerClass='SetOfAverages',
                      label="Input averages", important=True,
                      help='Select the input averages to be centered.')
        
        form.addParallelSection(threads=0, mpi=0)
    
    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.centerAveragesStep,
                                 self.inputAverages.get().getObjId())
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def centerAveragesStep(self, averagesId):
        inFn = self._getTmpPath('input_averages.mrcs')
        outFn = self._getStackFn()

        self.info("Writing input averages to: %s" % inFn)
        self.inputAverages.get().writeStack(inFn)

        self.runJob(self._getProgram('relion_image_handler'),
                    ' --shift_com --i %s --o %s' % (inFn, outFn))

    def createOutputStep(self):
        inputSet = self.inputAverages.get()
        avgSet = self._createSetOfAverages()

        avgSet.copyInfo(inputSet)
        avgSet.copyItems(inputSet, updateItemCallback=self._setFileName)
        self._defineOutputs(**{outputs.outputAverages.name: avgSet})
        self._defineTransformRelation(inputSet, avgSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        """ Just overwrite the default behaviour of the base class. """
        return []

    def _summary(self):
        summary = []

        if hasattr(self, "outputAverages"):
            summary.append('Class averages were aligned by relion_image_handler '
                           'using their center of mass.')

        return summary

    def _methods(self):
        return []

    def _getStackFn(self):
        return self._getPath('centered_averages.mrcs')

    def _setFileName(self, item, row):
        self._counter = getattr(self, '_counter', 0)
        self._counter += 1
        item.setLocation(self._counter, self._getStackFn())
