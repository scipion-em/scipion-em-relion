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

from enum import Enum

from pyworkflow.object import String, Integer
from pyworkflow.constants import PROD
from pyworkflow.protocol.params import (PointerParam, BooleanParam,
                                        IntParam, LabelParam)
from pwem.constants import ALIGN_PROJ
from pwem.protocols import ProtOperateParticles
from pwem.objects import SetOfParticles

import relion.convert as convert
from .protocol_base import ProtRelionBase


class outputs(Enum):
    outputParticles = SetOfParticles


class ProtRelionSubtract(ProtOperateParticles, ProtRelionBase):
    """
    Performs signal subtraction on cryo-EM particle images using RELION-based
    projection subtraction methods. The protocol removes selected density
    contributions from experimental particles in order to isolate flexible,
    heterogeneous, or weakly represented structural regions for downstream
    analysis.

    AI Generated:

    Signal Subtraction (ProtRelionSubtract) — User Manual
        Overview

        The Signal Subtraction protocol is designed to computationally remove
        unwanted or dominant structural features from cryo-EM particle images.
        This operation allows users to focus subsequent analyses on specific
        regions of interest, particularly flexible domains, secondary binding
        partners, or conformationally variable components that may otherwise be
        masked by stronger signal contributions.

        In practical cryo-EM workflows, signal subtraction is commonly used
        before focused classification, local refinement, or heterogeneity
        analysis. By subtracting projections of a known volume from the
        experimental particles, the remaining images emphasize the density that
        was intentionally preserved by the masking strategy. This often improves
        the interpretability of flexible regions and enables separation of
        biologically meaningful conformational states.

        Inputs and Biological Context

        The protocol can operate either from an existing RELION refinement or
        classification workflow or from externally provided particles and maps.
        In both cases, the particles must already contain projection alignment
        information because accurate orientation parameters are essential for
        generating projections that correctly match the experimental images.

        When using an existing RELION refinement, the protocol directly reuses
        the refinement geometry and associated metadata. This mode is typically
        preferred because it guarantees consistency between the reconstructed
        map and the particle orientations used during subtraction.

        Alternatively, users may provide an external reconstructed volume and a
        set of aligned particles. In this scenario, the biological validity of
        the subtraction strongly depends on the consistency between the map and
        the particle dataset. The reconstructed density should originate from
        the same particles and ideally from the same refinement strategy to
        preserve compatible intensity scaling and orientation conventions.

        Masking Strategy and Biological Interpretation

        The subtraction mask is the most biologically important component of
        the workflow because it defines which regions are retained and which
        are removed from the particles. Unlike masking strategies used for
        refinement, this protocol expects the mask to preserve the density of
        interest while excluding the signal intended for subtraction.

        In biological applications, masks are often designed to isolate mobile
        domains, ligand-binding regions, membrane-associated segments, or
        compositional variants within large assemblies. Careful mask design is
        essential because poor masking may introduce subtraction artifacts,
        residual density contamination, or unintended removal of meaningful
        structural features.

        Soft mask boundaries are generally recommended because abrupt edges may
        produce unrealistic Fourier artifacts that interfere with downstream
        classification or interpretation. For highly flexible systems, focusing
        the mask on structurally stable regions usually provides more reliable
        subtraction behavior.

        Particle Centering and Coordinate Handling

        The protocol supports optional recentering strategies after subtraction.
        Recentered particles can simplify downstream focused analyses by moving
        the retained region closer to the center of the particle box.

        One approach uses the center of mass of the masking region, which is
        often appropriate when the retained density forms a compact structural
        domain. Alternatively, users may provide explicit coordinates when prior
        biological knowledge suggests a more suitable reference position.

        Rewindowing into a smaller box size is also supported. This is useful
        when the retained signal occupies only a limited fraction of the
        original particle area, reducing computational cost during later
        refinement or classification stages. However, excessively aggressive
        cropping may truncate flexible regions or introduce edge artifacts.

        CTF Considerations

        The protocol includes optional handling of contrast transfer function
        effects during subtraction. Proper CTF treatment is important because
        projection subtraction occurs directly in the image domain and therefore
        depends on accurate modeling of microscope imaging effects.

        In most biological workflows, enabling CTF correction improves the
        quality of the subtraction and reduces residual artifacts. Additional
        options allow users to ignore low-frequency regions before the first
        CTF peak when the low-resolution CTF model is considered unreliable.
        Although potentially useful in difficult datasets, this strategy should
        generally be applied cautiously because it may alter the balance of
        structural signal across frequency ranges.

        Outputs and Downstream Applications

        The protocol produces a new particle dataset in which the selected
        structural signal has been computationally removed. The resulting
        particles preserve alignment information and can therefore be used
        directly in downstream focused classification, masked refinement,
        variability analysis, or local structural interpretation workflows.

        Biologically, the resulting particles often reveal heterogeneity that
        was previously obscured by dominant structural components. This is
        particularly valuable for studying flexible domains, transient
        interactions, partial ligand occupancy, or dynamic assemblies.

        In some cases, subtraction can also improve classification sensitivity
        by reducing the influence of large rigid cores that would otherwise
        dominate similarity measurements between particles.

        Practical Recommendations

        For most cryo-EM applications, the quality of the subtraction depends
        primarily on the accuracy of the alignment parameters and the biological
        relevance of the masking strategy. Before performing subtraction, users
        should ensure that the refinement producing the reference map is stable
        and well converged.

        It is generally advisable to begin with conservative masks and visually
        inspect the resulting particles for subtraction artifacts. Overly large
        masks or inaccurate centering strategies may remove meaningful signal
        or introduce distortions that complicate interpretation.

        Focused subtraction is particularly powerful when combined with focused
        classification without alignment, as this allows subtle structural
        variability to emerge without interference from dominant rigid regions.

        Final Perspective

        Signal subtraction is one of the most important focused-analysis tools
        in modern cryo-EM image processing because it enables direct study of
        structural heterogeneity within complex macromolecular assemblies.
        Successful application depends on biologically meaningful masking,
        reliable orientation information, and careful interpretation of the
        resulting particle populations. When properly applied, subtraction can
        reveal conformational variability and compositional differences that
        are otherwise inaccessible in conventional global refinements.
    """
    _label = 'subtract projection'
    _devStatus = PROD
    _possibleOutputs = outputs

    def _initialize(self):
        self._createFilenameTemplates()
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
                  'input_star': self._getExtraPath('input_particles.star'),
                  'output_star': self._getExtraPath('particles_subtracted.star')
                  }
        self._updateFilenamesDict(myDict)
    
    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('relionInput', BooleanParam,
                      default=True, important=True,
                      label="Start from Relion protocol?",
                      help="Set to Yes if you wish to use as input "
                           "a Relion protocol. Otherwise set it to No")
        form.addParam('inputProtocol', PointerParam,
                      important=True,
                      pointerClass='ProtRelionRefine3D, ProtRelionClassify3D,'
                                   'ProtRelionMultiBody',
                      label="Input Relion protocol",
                      condition="relionInput",
                      help="Select the 3D refinement/classification or multi-body "
                           "run which you want to use for subtraction. It will "
                           "use the maps from this run for the subtraction.")
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input map to be projected",
                      important=True,
                      condition="not relionInput",
                      help='Provide the input volume that will be used to '
                           'calculate projections, which will be subtracted '
                           'from the experimental particles. Make sure this '
                           'map was calculated by RELION from the same '
                           'particles as above, and preferably with those '
                           'orientations, as it is crucial that the absolute '
                           'greyscale is the same as in the experimental '
                           'particles.')
        form.addParam('useAll', BooleanParam, default=True,
                      label="Use all particles from input protocol?",
                      condition="relionInput",
                      help="If No, then you need to provide a subset of "
                           "particles below.")

        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      condition='not useAll',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles subset",
                      help='Select the particles which are a SUBSET of the '
                           'input protocol provided above.')

        form.addParam('inputParticlesAll', PointerParam,
                      pointerClass='SetOfParticles',
                      condition="not relionInput",
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Select the input particles.')

        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask of the signal to keep',
                      help="Provide a soft mask where the protein density "
                           "you wish to subtract from the experimental "
                           "particles is black (0) and the density you "
                           "wish to keep is white (1).\n"
                           "That is: *the mask should INCLUDE the part of the "
                           "volume that you wish to KEEP.*")

        form.addParam('saveFloat16', BooleanParam, default=True,
                      label="Write output in float16?",
                      help="Relion can write output images in float16 "
                           "MRC (mode 12) format to save disk space. "
                           "By default, float32 format is used.")

        form.addSection('Centering')
        form.addParam('help1', LabelParam,
                      condition="not relionInput",
                      label="This section is only used if "
                            "starting from Relion input.")
        form.addParam('centerOnMask', BooleanParam, default=True,
                      label="Do center subtracted images on mask?",
                      condition="relionInput",
                      help="If set to Yes, the subtracted particles will "
                           "be centered on projections of the "
                           "center-of-mass of the input mask.")
        form.addParam('centerOnCoord', BooleanParam, default=False,
                      condition='(not centerOnMask) and relionInput',
                      label="Do center on my coordinates?",
                      help="If set to Yes, the subtracted particles will "
                           "be centered on projections of the x,y,z "
                           "coordinates below. The unit is pixel, not "
                           "angstrom. The origin is at the center of the box, "
                           "not at the corner.")

        line = form.addLine('Center coordinate (px)',
                            condition='centerOnCoord and relionInput',
                            help='Coordinate of the 3D center (in pixels).')
        line.addParam('cX', IntParam, default=0, condition='centerOnCoord and relionInput',
                      label='X')
        line.addParam('cY', IntParam, default=0, condition='centerOnCoord and relionInput',
                      label='Y')
        line.addParam('cZ', IntParam, default=0, condition='centerOnCoord and relionInput',
                      label='Z')

        form.addParam('newBoxSize', IntParam, default=-1,
                      condition="relionInput",
                      label="New box size",
                      help="Provide a non-negative value to re-window the "
                           "subtracted particles in a smaller box size.")

        form.addSection(label='CTF')
        form.addParam('help', LabelParam,
                      condition="relionInput",
                      label="This section is only used if "
                            "NOT starting from Relion protocol.")
        form.addParam('doCTF', BooleanParam, default=True,
                      label='Do CTF-correction?',
                      condition="not relionInput",
                      help='If set to Yes, CTFs will be corrected inside the '
                           'MAP refinement. The resulting algorithm '
                           'intrinsically implements the optimal linear, or '
                           'Wiener filter. Note that input particles should '
                           'contains CTF parameters.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam, default=False,
                      label='Ignore CTFs until first peak?',
                      condition="not relionInput",
                      help='If set to Yes, then CTF-amplitude correction will '
                           'only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF '
                           'model is inadequate at the lowest resolution. '
                           'Still, in general using higher amplitude contrast '
                           'on the CTFs (e.g. 10-20%) often yields better '
                           'results. Therefore, this option is not generally '
                           'recommended.')

        form.addParallelSection(threads=0, mpi=2)
    
    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self.isRelionInput = self.relionInput.get()
        self._initialize()

        if not self.useAll or not self.isRelionInput:
            self._insertFunctionStep(self.convertInputStep, needsGPU=False)

        self._insertFunctionStep(self.subtractStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)
    
    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        """ Write the input images as a Relion star file. """
        if self.isRelionInput:
            imgSet = self.inputParticles.get()
        else:
            imgSet = self.inputParticlesAll.get()

        convert.writeSetOfParticles(
            imgSet, self._getFileName('input_star'),
            outputDir=self._getExtraPath(), alignType=ALIGN_PROJ)

    def subtractStep(self):
        if self.isRelionInput:
            self.subtractStepRelion()
        else:
            self.subtractStepNoRelion()

    def subtractStepNoRelion(self):
        volume = self.inputVolume.get()
        volFn = convert.convertBinaryVol(volume,
                                         self._getExtraPath())
        params = ' --i %s --subtract_exp' % volFn
        params += ' --angpix %0.3f' % volume.getSamplingRate()
        params += self._convertMask(resize=False, invert=True)

        if self.doCTF:
            params += ' --ctf'
            if self.ignoreCTFUntilFirstPeak:
                params += ' --ctf_intact_first_peak'
            if self._getInputParticles().isPhaseFlipped():
                params += ' --ctf_phase_flip'

        params += ' --ang %s  --o %s' % (
            self._getFileName('input_star'),
            self._getFileName('output_star').replace(".star", ""))

        self.runJob('relion_project', params)

    def subtractStepRelion(self):
        inputProt = self.inputProtocol.get()
        inputProt._initialize()

        fnOptimiser = inputProt._getOptimiserFile()
        params = " --i %s --o %s --new_box %s" % (fnOptimiser,
                                                  self._getExtraPath(),
                                                  self.newBoxSize.get())

        if not self.useAll:
            params += " --data %s" % self._getFileName('input_star')
        if self.centerOnMask:
            params += " --recenter_on_mask"
        elif self.centerOnCoord:
            params += " --center_x %d --center_y %d --center_z %d" % (
                self.cX, self.cY, self.cZ)

        if self.saveFloat16:
            params += " --float16"

        params += self._convertMask()
        self.runJob(self._getProgram('relion_particle_subtract'), params)

    def createOutputStep(self):
        imgSet = self._getInputParticles()
        outImgSet = self._createSetOfParticles()
        outImgsFn = self._getFileName('output_star')
        outImgSet.copyInfo(imgSet)
        outImgSet.setAlignmentProj()

        px = imgSet.getSamplingRate()
        self.reader = convert.createReader(alignType=ALIGN_PROJ,
                                           pixelSize=px)
        mdIter = convert.Table.iterRows('particles@' + outImgsFn,
                                        types=convert.LABELS_DICT)
        outImgSet.copyItems(imgSet, doClone=False,
                            updateItemCallback=self._updateItem,
                            itemDataIterator=mdIter)

        self._defineOutputs(**{outputs.outputParticles.name: outImgSet})
        self._defineTransformRelation(imgSet, outImgSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if not self.useAll:
            self._validateDim(self.inputParticles(),
                              self._getInputParticles().getXDim(),
                              errors, 'Input particles subset',
                              'Input particles from 3D protocol')
        if self.numberOfMpi > 1 and (not self.relionInput.get()):
            errors.append("Use of several CPUs when input is not relion "
                          "protocol is not supported")

        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output is not ready yet.")
        else:
            summary.append('Projections of the masked input volume were '
                           'subtracted from original particles.')

        return summary
    
    # -------------------------- UTILS functions ------------------------------
    def _updateItem(self, particle, row):
        if self.isRelionInput:
            # FIXME: check if other attrs need saving
            particle._rlnRandomSubset = Integer(row.rlnRandomSubset)
            self.reader.setParticleTransform(particle, row)
        particle._rlnImageOriginalName = String(row.rlnImageOriginalName)
        newFn = row.rlnImageName
        newLoc = convert.relionToLocation(newFn)
        particle.setLocation(newLoc)

    def _getInputParticles(self):
        if self.isRelionInput:
            inputProt = self.inputProtocol.get()
            return inputProt.outputParticles
        else:
            return self.inputParticlesAll.get()

    def _convertMask(self, invert=False, resize=True):
        tmp = self._getTmpPath()
        if resize:
            newDim = self._getInputParticles().getXDim()
            newPix = self._getInputParticles().getSamplingRate()
        else:
            newDim = None
            newPix = None
        maskFn = convert.convertMask(self.refMask.get(),
                                     tmp, newPix, newDim, invert=invert)
        return ' --mask %s' % maskFn
