# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
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

import pyworkflow.protocol.params as params
from pwem.protocols import ProtAnalysis3D

from .protocol_base import ProtRelionBase


class ProtRelionPlotBfactor(ProtAnalysis3D, ProtRelionBase):
    """
    Relion reconstruction protocol to produce ResLog-like plot
    """
    _label = 'plot b-factor'
    IS_3D = True

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'finalVolume': self._getExtraPath("relion_class001.mrc"),
        }
        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('protRefine', params.PointerParam,
                      important=True,
                      pointerClass="ProtRefine3D",
                      label='Input refinement protocol',
                      help='Select any previous refinement protocol '
                           'with all particles. ')
        form.addParam('protPostProcess', params.PointerParam,
                      important=True,
                      pointerClass="ProtRelionPostprocess",
                      label='Input postprocess protocol',
                      help='Select any previous postprocess protocol, '
                           'which is required for resolution assessment.')

        line = form.addLine('Number of particles: ',
                            expertLevel=params.LEVEL_ADVANCED,
                            help='Specify minimum and maximum number of '
                                 'particles.')
        line.addParam('minPtcls', params.IntParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default='100', label='min')
        line.addParam('maxPtcls', params.IntParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default='9999999', label='max')

        form.addParam('initialLowPassFilterA', params.FloatParam,
                      default=60,
                      label='Initial low-pass filter (A)',
                      help='It is recommended to strongly low-pass filter your '
                           'initial reference map. If it has not yet been '
                           'low-pass filtered, it may be done internally using '
                           'this option. If set to 0, no low-pass filter will '
                           'be applied to the initial reference(s).')

        form.addSection('Compute')
        self._defineComputeParams(form)

        form.addParallelSection(threads=1, mpi=3)

    # -------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        objId = self.protRefine.get().getObjId()
        self._createFilenameTemplates()
        self._defineParamDict()
        self._insertFunctionStep('convertInputStep', objId)
        self._insertFunctionStep('runRelionStep', self.paramDict)
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, protId):
        pass

    def runRelionStep(self, paramDict):
        '''
        0. Get rlnFinalResolution and rlnBfactorUsedForSharpening from input postprocess job for all ptcls.
        1. Create subset of input ptcls of size MIN, randomly
        2. Run refine 3d.
        3. Run postprocess. Get resolution and b-factor.
        4. Double number of ptcls. Repeat 1-3 while N<MAX and N<total
        5. Make ResLog + b-factor plots.
        '''
        params = ' '.join(['%s %s' % (k, str(v))
                           for k, v in self.paramDict.items()])

        program = 'relion_postprocess'
        if self.numberOfMpi > 1:
            program += '_mpi'

        self.runJob(program, params)

    def createOutputStep(self):
        pass

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        errors = []

        return errors

    def _citations(self):
        return []

    def _summary(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        summary = []

        return summary

    # -------------------------- UTILS functions ------------------------------
    def _defineParamDict(self):
        """ Define all parameters to run relion_refine"""

        self.paramDict = {}
