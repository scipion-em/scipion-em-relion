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

import os
from emtable import Table

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.protocols import ProtAnalysis3D
from pwem.objects import FSC, SetOfFSCs
from pwem.emlib.image import ImageHandler

from ..constants import (FSC_TYPE_OVERALL, FSC_TYPE_MODEL_MAP,
                         FSC_TYPE_WORK_FREE)
from .protocol_base import ProtRelionBase


class ProtRelionCalculateFSC(ProtAnalysis3D, ProtRelionBase):
    """
    Relion protocol to calculate various FSC curves using relion_image_handler.
    """
    _label = 'calculate fsc'
    _devStatus = PROD
    _possibleOutputs = {
        'outputFSC': FSC,
        'outputSetOfFSCs': SetOfFSCs
    }

    def _createFilenameTemplates(self):
        """ Centralize how the files are called. """
        myDict = {
            'model': self._getTmpPath("input_model_final.mrc"),
            'model_half1': self._getTmpPath("input_model_half1.mrc"),
            'half1': self._getTmpPath("input_half1_unfil.mrc"),
            'half2': self._getTmpPath("input_half2_unfil.mrc"),
            'map': self._getTmpPath("input_map.mrc"),
            # outputs
            'fsc_overall': self._getExtraPath("fsc_overall.star"),
            'fsc_model-map': self._getExtraPath("fsc_model-map.star"),
            'fsc_work': self._getExtraPath("fsc_work.star"),
            'fsc_free': self._getExtraPath("fsc_free.star")
        }
        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('fscType', params.EnumParam, default=0,
                      choices=['FSC overall', 'FSC model-map',
                               'FSC work / FSC free'],
                      label="Select FSC type to compute",
                      help="1) FSC overall - between two half-maps\n"
                           "2) FSC model-map - between atomic model and refined map\n"
                           "3) FSC work - model refined against half-map 1, "
                           "compared to half-map 1\n"
                           "FSC free - model refined against half-map 1, "
                           "compared to half-map 2")

        form.addParam('half1', params.PointerParam, pointerClass='Volume',
                      condition="fscType!=%d" % FSC_TYPE_MODEL_MAP,
                      label="Input half map 1",
                      important=True, allowsNull=True)
        form.addParam('half2', params.PointerParam, pointerClass='Volume',
                      condition="fscType!=%d" % FSC_TYPE_MODEL_MAP,
                      label="Input half map 2",
                      important=True, allowsNull=True)

        form.addParam('map', params.PointerParam, pointerClass='Volume',
                      condition="fscType==%d" % FSC_TYPE_MODEL_MAP,
                      label="Final map",
                      important=True, allowsNull=True)

        form.addParam('model', params.PointerParam, pointerClass='AtomStruct',
                      condition="fscType==%d" % FSC_TYPE_MODEL_MAP,
                      label="Final atomic model",
                      important=True, allowsNull=True)

        form.addParam('model_half1', params.PointerParam, pointerClass='AtomStruct',
                      condition="fscType==%d" % FSC_TYPE_WORK_FREE,
                      label="Atomic model refined against half-map 1",
                      important=True, allowsNull=True)

        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('calculateFSCStep')
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, *args):
        """ Create links and convert pdb to mrc. """
        if self._getFSCType() == FSC_TYPE_OVERALL:
            self._createInputLink(self.half1.get().getFileName(), 'half1')
            self._createInputLink(self.half2.get().getFileName(), 'half2')

        elif self._getFSCType() == FSC_TYPE_MODEL_MAP:
            self._createInputLink(self.map.get().getFileName(), 'map')

            box = self.map.get().getDim()[0]
            apix = self.map.get().getSamplingRate()

            self._model2Map(self.model.get().getFileName(),
                            self._getFileName('model'), box, apix)

        else:
            self._createInputLink(self.half1.get().getFileName(), 'half1')
            self._createInputLink(self.half2.get().getFileName(), 'half2')

            box = self.half1.get().getDim()[0]
            apix = self.half1.get().getSamplingRate()

            self._model2Map(self.model_half1.get().getFileName(),
                            self._getFileName('model_half1'), box, apix)

    def calculateFSCStep(self):
        """ Run relion_image_handler. """
        if self._getFSCType() == FSC_TYPE_OVERALL:
            angpix = self.half1.get().getSamplingRate()
            params = self._getParams(self._getFileName('half1'),
                                     self._getFileName('half2'),
                                     self._getFileName('fsc_overall'), angpix)
        elif self._getFSCType() == FSC_TYPE_MODEL_MAP:
            angpix = self.map.get().getSamplingRate()
            params = self._getParams(self._getFileName('model'),
                                     self._getFileName('map'),
                                     self._getFileName('fsc_model-map'), angpix)
        else:
            # FSC work
            angpix = self.half1.get().getSamplingRate()
            params = self._getParams(self._getFileName('model_half1'),
                                     self._getFileName('half1'),
                                     self._getFileName('fsc_work'), angpix)
            self._runProgram('relion_image_handler', params)

            # FSC free
            params = self._getParams(self._getFileName('model_half1'),
                                     self._getFileName('half2'),
                                     self._getFileName('fsc_free'), angpix)

        self._runProgram('relion_image_handler', params)

    def createOutputStep(self):
        if self._getFSCType() == FSC_TYPE_OVERALL:
            fsc = self._getFSC(self._getFileName('fsc_overall'),
                               'FSC overall')
            self._defineOutputs(outputFSC=fsc)

        elif self._getFSCType() == FSC_TYPE_MODEL_MAP:
            fsc = self._getFSC(self._getFileName('fsc_model-map'),
                               'FSC model-map')
            self._defineOutputs(outputFSC=fsc)

        else:
            fscSet = self._createSetOfFSCs()
            fscw = self._getFSC(self._getFileName('fsc_work'),
                                'FSC work')
            fscf = self._getFSC(self._getFileName('fsc_free'),
                                'FSC free')
            fscSet.append(fscw)
            fscSet.append(fscf)
            self._defineOutputs(outputSetOfFSCs=fscSet)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []

        from pwem import Domain
        try:
            _ = Domain.importFromPlugin('eman2', doRaise=True)
        except:
            errors.append("EMAN2 is required to convert pdb to mrc map")
        return errors

    def _citations(self):
        return ['Amunts2014']

    def _summary(self):
        summary = []

        if hasattr(self, "outputFSC") or hasattr(self, "outputSetOfFSCs"):
            summary.append("FSC calculation ran with relion_image_handler")

        return summary

    # -------------------------- UTILS functions ------------------------------
    def _createInputLink(self, fn, link):
        if pwutils.getExt(fn) == ".mrc":
            pwutils.createAbsLink(os.path.abspath(fn), self._getFileName(link))
        else:
            ih = ImageHandler()
            ih.convert(fn, self._getFileName(link))

    def _getFSCType(self):
        return self.fscType.get()

    def _model2Map(self, model, map, box, apix):
        """ Convert model pdb to mrc map. """
        from pwem import Domain
        eman2 = Domain.importFromPlugin('eman2',
                                        errorMsg='EMAN2 is required to convert pdb to mrc map',
                                        doRaise=True)
        from pyworkflow.utils.process import runJob

        args = "%s %s --box %d --apix %0.3f" % (model, map, box, apix)
        runJob(self._log, eman2.Plugin.getProgram('e2pdb2mrc.py'), args,
               env=eman2.Plugin.getEnviron())

    def _getFSC(self, fn, label):
        fsc = FSC(objLabel=label)
        table = Table(fileName=fn, tableName='fsc')
        resolution_inv = table.getColumnValues('rlnResolution')
        frc = table.getColumnValues('rlnFourierShellCorrelation')
        fsc.setData(resolution_inv, frc)
        return fsc

    def _getParams(self, fn1, fn2, output, angpix):
        params = ' '.join([
            '--i %s ' % fn1,
            '--fsc %s' % fn2,
            '--angpix %0.3f' % angpix,
            '> %s' % output
        ])
        return params
