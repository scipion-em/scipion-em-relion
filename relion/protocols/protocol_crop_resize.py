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

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD
from pwem.protocols import ProtPreprocessVolumes
from pwem.emlib.image import ImageHandler
from pwem.objects import Volume


class ProtRelionResizeVolume(ProtPreprocessVolumes):
    """ This protocol rescales/resizes 3D volumes using relion_image_handler. """

    _label = 'crop/resize volumes'
    _devStatus = PROD

    def __init__(self, **kwargs):
        ProtPreprocessVolumes.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize how the protocol files are called. """
        myDict = {
            'input_vol': self._getTmpPath('volume_%(volId)03d.mrc'),
            'output_vol': self._getExtraPath('volume_rescaled_%(volId)03d.mrc')
        }
        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=pwutils.Message.LABEL_INPUT)
        form.addParam('inputVolumes', params.PointerParam, important=True,
                      label=pwutils.Message.LABEL_INPUT_VOLS,
                      pointerClass='Volume, SetOfVolumes',
                      help='Can be a Volume or a SetOfVolumes.')
        form.addParam('doRescale', params.BooleanParam, default=False,
                      label='Rescale volumes?')
        form.addParam('rescaleSamplingRate', params.FloatParam,
                      default=1.0,
                      condition='doRescale',
                      label='New sampling rate (Å/px)')
        form.addParam('doResize', params.BooleanParam, default=False,
                      label='Resize volumes to a new box?')
        form.addParam('resizeSize', params.IntParam, default=0,
                      condition='doResize',
                      label='New box size (px)',
                      help='Provide even box size.')

        form.addParallelSection(threads=0, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('resizeVolumesStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        vols = self.inputVolumes.get()
        self.convertedVols = []

        if isinstance(vols, Volume):
            fn = self._convertVol(vols, 1)
            self.convertedVols.append(fn)
        else:
            for i, vol in enumerate(vols):
                fn = self._convertVol(vol, i)
                self.convertedVols.append(fn)

    def resizeVolumesStep(self):
        argDict = {' --angpix': self.inputVolumes.get().getSamplingRate()}

        if self.doRescale:
            argDict[' --rescale_angpix'] = self.rescaleSamplingRate.get()

        if self.doResize:
            argDict[' --new_box'] = self.resizeSize.get()

        for i, fn in enumerate(self.convertedVols):
            argDict[' --i'] = fn
            argDict[' --o'] = self._getFileName('output_vol', volId=i+1)
            args = ' '.join(['%s %s' % (k, v) for k, v in argDict.items()])
            self.runJob("relion_image_handler", "".join(args))

    def createOutputStep(self):
        volInput = self.inputVolumes.get()
        if isinstance(volInput, Volume):
            # Create the output with the same class as the input, that should
            # be Volume or a subclass of Volume e.g. VolumeMask
            volClass = volInput.getClass()
            vol = volClass()
            vol.copyInfo(volInput)
            vol.setLocation(self._getFileName('output_vol', volId=1))
            vol.setSamplingRate(self._getNewSampling())
            self._defineOutputs(outputVol=vol)
        else:
            volumes = self._createSetOfVolumes()
            volumes.copyInfo(volInput)
            volumes.setSamplingRate(self._getNewSampling())
            for i, obj in enumerate(volInput.iterItems()):
                vol = obj
                vol.setLocation(self._getFileName('output_vol', volId=i+1))
                vol.setSamplingRate(self._getNewSampling())
                volumes.append(vol)
            self._defineOutputs(outputVol=volumes)

        self._defineTransformRelation(self.inputVolumes, self.outputVol)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if not self.doRescale and not self.doResize:
            errors.append("You have to select at least one option!")

        if self.doResize:
            if not self.resizeSize.get() % 2 == 0:
                errors.append("Only even box sizes are allowed!")

        return errors

    def _summary(self):
        messages = []

        if self.doRescale:
            messages.append("- Rescaled volumes to pixel size of %0.3f"
                            % self.rescaleSamplingRate.get())
        if self.doResize:
            messages.append("- Resized volumes to box size of %d"
                            % self.resizeSize.get())

        return messages

    # -------------------------- UTILS functions ------------------------------
    def _convertVol(self, vol, index):
        ih = ImageHandler()
        fn = vol.getFileName()
        if not fn.endswith('.mrc'):
            newFn = self._getFileName('input_vol', volId=index)
            ih.convert(fn, newFn)
            return newFn
        return fn

    def _getNewSampling(self):
        if self.doRescale:
            return self.rescaleSamplingRate.get()
        else:
            return self.inputVolumes.get().getSamplingRate()
