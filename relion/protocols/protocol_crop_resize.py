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
from pwem.convert import Ccp4Header


class ProtRelionResizeVolume(ProtPreprocessVolumes):
    """
    Rescales and resizes 3D cryo-EM volumes while preserving their
    structural content and compatibility with downstream image analysis
    workflows. The protocol is designed to adapt map dimensions and voxel
    sizes to match the requirements of reconstruction, refinement,
    visualization, or comparative structural analysis.

    AI Generated:

    Resize Volume (ProtRelionResizeVolume) — User Manual
        Overview

        The Resize Volume protocol modifies the sampling rate and box size
        of three-dimensional cryo-EM maps. Its primary purpose is to make
        volumes compatible with downstream processing steps, facilitate
        comparisons between datasets, or optimize computational efficiency.
        In cryo-EM workflows, resizing and rescaling are common operations
        when combining maps from different sources, preparing references for
        refinement, or adapting reconstructions to specific software
        requirements.

        For biological users, this protocol is particularly useful when maps
        have incompatible voxel sizes or box dimensions. Standardizing these
        properties allows reconstructions to be compared directly, aligned
        consistently, or processed together in later stages of analysis.
        Although resizing is mathematically straightforward, it can strongly
        influence map interpretation and should therefore be applied with
        biological caution.

        Inputs and General Workflow

        The protocol accepts either a single volume or an entire set of
        volumes. This flexibility is useful in practical cryo-EM projects,
        where users may need to process individual reconstructions, multiple
        conformational states, masks, or collections of maps originating from
        classification procedures.

        The workflow may involve two related but distinct operations:
        rescaling the voxel size and resizing the box dimensions. These
        operations can be performed independently or together, depending on
        the biological and computational objective of the analysis.

        Rescaling Volumes

        Rescaling changes the effective sampling rate of the map. This
        operation is commonly required when volumes from different datasets
        or software packages must be brought into a consistent physical
        scale. For example, a reconstruction generated with one pixel size
        may need to match the sampling of another reconstruction before
        comparison, subtraction, docking, or flexible fitting.

        From a biological perspective, preserving consistent spatial scaling
        is critical because structural measurements, molecular dimensions,
        and atomic interpretation all depend directly on the voxel size.
        Incorrect scaling can lead to inaccurate structural conclusions or
        failed downstream refinements.

        Users should also recognize that rescaling may introduce small
        adjustments due to numerical constraints associated with even box
        dimensions. As a result, the final voxel size may differ slightly
        from the requested value. These differences are usually minor but
        should still be considered when performing quantitative structural
        analysis.

        Resizing the Box

        Resizing changes the number of voxels defining the three-dimensional
        map volume. This operation is often used to crop unnecessary solvent
        regions, reduce computational cost, or adapt maps to software
        requirements that expect specific box dimensions.

        Smaller boxes reduce memory usage and accelerate downstream
        calculations, which is especially valuable during iterative
        refinement or large-scale processing. However, overly aggressive
        cropping can remove biologically important density, truncate flexible
        regions, or introduce edge artifacts that negatively affect later
        analysis.

        Enlarging the box may also be useful in some workflows, particularly
        when additional surrounding space is needed for alignment,
        visualization, or flexible motion analysis. In all cases, users
        should ensure that the molecular complex remains fully contained
        within the resized volume.

        Handling Half-Maps

        The protocol supports processing of half-maps together with the main
        reconstruction. Maintaining consistency between consensus maps and
        corresponding half-maps is biologically important because many
        validation procedures, including resolution estimation and local
        refinement assessment, rely on the correct relationship between
        these maps.

        By preserving compatible dimensions and sampling across all related
        reconstructions, the protocol helps maintain reliable downstream
        validation and interpretation.

        Outputs and Their Interpretation

        After execution, the protocol produces resized or rescaled volumes
        with updated voxel sizes and box dimensions. The resulting outputs
        retain the structural information of the original maps while being
        adapted to the desired computational or biological context.

        For sets of volumes, each output preserves its original identity and
        metadata while sharing the newly standardized dimensions and
        sampling. This is especially useful for comparative structural
        studies involving multiple conformational states or reconstruction
        methods.

        Practical Recommendations

        In routine cryo-EM workflows, users should carefully evaluate
        whether resizing or rescaling is biologically justified before
        applying the protocol. Rescaling is appropriate when maps originate
        from different pixel sizes, whereas resizing is mainly useful for
        computational optimization or compatibility.

        When cropping volumes, it is advisable to preserve sufficient solvent
        margins around the molecular complex in order to avoid truncation
        artifacts. Flexible domains, membrane regions, or extended filament
        segments may require larger boxes than initially expected.

        For high-resolution work, users should verify that interpolation
        effects introduced during rescaling do not compromise fine structural
        details. Visual inspection of the resized maps is strongly
        recommended before continuing with refinement or interpretation.

        Final Perspective

        Volume resizing and rescaling are fundamental preparation steps in
        cryo-EM image analysis. Although these operations are computational
        in nature, they have direct biological consequences because they
        determine how structural information is represented, compared, and
        interpreted throughout downstream workflows. Careful control of voxel
        size and box dimensions helps ensure reliable structural analysis and
        robust integration of cryo-EM datasets.
    """

    _label = 'crop/resize volumes'
    _devStatus = PROD

    def __init__(self, **kwargs):
        ProtPreprocessVolumes.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize how the protocol files are called. """
        myDict = {
            'input_vol': self._getTmpPath('volume_%(volId)03d.mrc'),
            'output_vol': self._getExtraPath('volume_rescaled_%(volId)03d.mrc'),
            'output_half1': self._getExtraPath('half1_rescaled_%(volId)03d.mrc'),
            'output_half2': self._getExtraPath('half2_rescaled_%(volId)03d.mrc')
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
                      label='New sampling rate (A/px)')
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
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.resizeVolumesStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

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
        inputData = self.inputVolumes.get()
        argDict = {' --angpix': inputData.getSamplingRate()}

        if self.doRescale:
            argDict[' --rescale_angpix'] = self.rescaleSamplingRate.get()

        if self.doResize:
            argDict[' --new_box'] = self.resizeSize.get()

        if isinstance(inputData, Volume):
            self.runResizeCmd(self.convertedVols, argDict)
            if inputData.hasHalfMaps():
                halves = inputData.getHalfMaps().split(',')
                self.convertedVols = []
                self.convertedVols.append(halves[0])
                self.runResizeCmd(self.convertedVols, argDict, keyVol='output_half1')

                self.convertedVols = []
                self.convertedVols.append(halves[1])
                self.runResizeCmd(self.convertedVols, argDict, keyVol='output_half2')
        else:
            self.runResizeCmd(self.convertedVols, argDict)

    def runResizeCmd(self, convertedVols, argDict, keyVol='output_vol'):
        for i, fn in enumerate(convertedVols):
            argDict[' --i'] = fn
            argDict[' --o'] = self._getFileName(keyVol, volId=i+1)
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
            if volInput.hasHalfMaps():
                halves = [self._getFileName('output_half1', volId=1), self._getFileName('output_half2', volId=1)]
                vol.setHalfMaps(halves)
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

        if self.doResize and not (self.resizeSize.get() % 2 == 0):
            errors.append("Only even box sizes are allowed!")

        return errors

    def _summary(self):
        messages = []

        if hasattr(self, "outputVol"):
            if self.doResize:
                messages.append(f"Resized volumes to box size of "
                                f"*{self.resizeSize.get()}* px")
            if self.doRescale:
                messages.append(f"The output pixel size might be different than the "
                                f"requested {self.rescaleSamplingRate.get()}A "
                                f"due to rounding of the box size "
                                f"to an even number by RELION")
        else:
            messages.append("Output is not ready")

        return messages

    # -------------------------- UTILS functions ------------------------------
    def _convertVol(self, vol, index):
        fn = vol.getFileName()
        self._convertFnVol(fn, index)
        return fn

    def _convertFnVol(self, fn, index):
        ih = ImageHandler()
        if not fn.endswith('.mrc'):
            newFn = self._getFileName('input_vol', volId=index)
            ih.convert(fn, newFn)
            return newFn

    def _getNewSampling(self):
        if self.doRescale:
            ccp4header = Ccp4Header(self._getFileName('output_vol', volId=1),
                                    readHeader=True)
            sampling, _, _ = ccp4header.getSampling()

            return sampling
        else:
            return self.inputVolumes.get().getSamplingRate()

