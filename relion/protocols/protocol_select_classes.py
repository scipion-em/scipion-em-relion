# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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
from emtable import Table

from pyworkflow.object import Float
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.protocols import ProtProcessParticles
from pwem.objects import SetOfClasses2D, SetOfParticles

from relion import Plugin
from .protocol_base import ProtRelionBase
from relion.convert import locationToRelion


class outputs(Enum):
    outputClasses = SetOfClasses2D
    outputParticles = SetOfParticles


class ProtRelionSelectClasses2D(ProtProcessParticles, ProtRelionBase):
    """
    Relion protocol to auto-select 2D class averages.
    """
    _label = '2D class ranker'
    _devStatus = PROD
    _possibleOutputs = outputs

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputProtocol', params.PointerParam,
                      pointerClass='ProtRelionClassify2D',
                      label="Input Relion 2D classification",
                      important=True)
        form.addParam('minThreshold', params.FloatParam, default=0.5,
                      label='Min. threshold for auto-selection',
                      help='Only classes with a predicted threshold '
                           'above this value will be selected.')
        form.addParam('minParts', params.IntParam, default=-1,
                      label='Select at least this many particles',
                      help='Even if they have scores below the minimum '
                           'threshold, select at least this many particles '
                           'with the best scores.')
        form.addParam('minCls', params.IntParam, default=-1,
                      label='OR: Select at least this many classes',
                      help='Even if they have scores below the minimum '
                           'threshold, select at least this many classes '
                           'with the best scores.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSelectStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def runSelectStep(self):
        inputProt = self.inputProtocol.get()
        inputProt._initialize()
        fnOptimiser = inputProt._getOptimiserFile()
        params = " --opt %s --o %s --min_score %s" % (fnOptimiser,
                                                      self._getExtraPath(),
                                                      self.minThreshold.get())
        params += " --fn_sel_parts particles.star"
        params += " --fn_sel_classavgs class_averages.star"
        params += " --fn_root rank --do_granularity_features"
        params += " --auto_select"

        if not Plugin.IS_GT50():
            params += " --python $CONDA_PREFIX/bin/python"

        if self.minParts != -1:
            params += " --select_min_nr_particles %d" % self.minParts
        if self.minCls != -1:
            params += " --select_min_nr_classes %d" % self.minCls

        self.runJob("%s && relion_class_ranker" % Plugin.getActivationCmd(), params)

    def createOutputStep(self):
        table = Table(fileName=self._getExtraPath('backup_selection.star'))
        selected = len([s for s in table.getColumnValues('rlnSelected') if s])

        if selected:
            classesStar = self._getExtraPath('class_averages.star')
            clsDict = {row.rlnReferenceImage: row
                       for row in Table.iterRows(classesStar)}

            inputClasses = self.inputProtocol.get().outputClasses
            outputClasses = SetOfClasses2D.create(self._getExtraPath())
            outputClasses.copyInfo(inputClasses)
            inputParticles = inputClasses.getImages()
            outputParticles = SetOfParticles.create(self._getExtraPath())
            outputParticles.copyInfo(inputParticles)

            def _getClassRow(cls2d):
                idx, fn = cls2d.getRepresentative().getLocation()
                return clsDict.get(locationToRelion(idx, fn), None)

            def _updateClass(cls2d):
                row = _getClassRow(cls2d)
                cls2d._rlnPredictedClassScore = Float(row.rlnPredictedClassScore)
                cls2d._rlnEstimatedResolution = Float(row.rlnEstimatedResolution)

            outputClasses.appendFromClasses(inputClasses,
                                            filterClassFunc=_getClassRow,
                                            updateClassCallback=_updateClass)

            self.summaryVar.set(f"Selected *{selected}* best classes.\n"
                                f"Threshold: *{self.minThreshold.get()}*")
            outputParticles.appendFromClasses(outputClasses)
            self._defineOutputs(**{outputs.outputClasses.name: outputClasses,
                                   outputs.outputParticles.name: outputParticles})
            self._defineSourceRelation(inputClasses, outputClasses)
            self._defineSourceRelation(inputParticles, outputParticles)
        else:
            self.summaryVar.set("No classes were selected.\n"
                                "Try with a lower threshold.")

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        return [self.summaryVar.get(default="No summary information")]

    def _validate(self):
        errors = []
        if self.minParts != -1 and self.minCls != -1:
            errors.append("You cannot choose both min. number of particles "
                          "and classes.")

        return errors
