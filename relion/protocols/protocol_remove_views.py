# ******************************************************************************
# *
# * Authors:     Grigory Sharov     (gsharov@mrc-lmb.cam.ac.uk)
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
# ******************************************************************************

import numpy as np
from random import sample
from enum import Enum

import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.protocols import ProtParticles
from pwem.objects import SetOfParticles
import pwem.convert.transformations as tfs


class outputs(Enum):
    outputParticles = SetOfParticles


class ProtRelionRemovePrefViews(ProtParticles):
    """
    Removes preferential particle orientations from a cryo-EM particle set in
    order to reduce angular bias and improve the isotropy of downstream
    reconstructions. The protocol enables selective elimination of particles
    whose viewing directions fall within user-defined angular regions of the
    Euler sphere. It is particularly useful when datasets contain excessive
    representation of specific orientations that may negatively affect 3D
    refinement, directional resolution, or map interpretability.

    Inspired by https://github.com/leschzinerlab/Relion

    AI Generated:

    Remove Preferential Views (ProtRelionRemovePrefViews) - User Manual
        Overview

        The Remove Preferential Views protocol is designed to reduce angular
        overrepresentation in cryo-EM particle datasets by selectively removing
        particles that belong to preferred orientations. In many single-particle
        cryo-EM experiments, molecules adsorb to the air-water interface or
        support film in a non-random manner, producing datasets dominated by a
        limited set of views. Although large numbers of particles may initially
        appear beneficial, severe angular imbalance can reduce directional
        resolution, introduce anisotropy, and compromise the quality of the
        reconstructed map.

        This protocol helps users create a more balanced angular distribution by
        removing particles located within specified orientation ranges. The
        resulting dataset is often better suited for testing reconstruction
        robustness, evaluating angular coverage, or preparing datasets for more
        isotropic refinement strategies.

        Biological Motivation

        Preferential orientation is one of the most common limitations in
        cryo-EM data collection. Biological assemblies with flat surfaces,
        membrane-associated regions, or strong adsorption tendencies frequently
        adopt only a narrow subset of orientations on the grid. As a result,
        certain projection directions become heavily oversampled while others
        remain poorly represented or completely absent.

        Removing part of the dominant orientations can sometimes improve the
        effective balance of the dataset and reveal whether reconstruction
        artifacts are linked to angular bias. This strategy is particularly
        useful during method development, refinement validation, or when testing
        the contribution of underrepresented orientations to the final map.

        Inputs and Orientation Requirements

        The protocol requires a particle set containing projection alignment
        information. Since particle orientations are defined in angular space,
        the input dataset must already include valid rotational assignments
        obtained from previous refinement or classification procedures.

        The protocol operates on Euler angles commonly described by rot, tilt,
        and optionally psi. Rot and tilt define the viewing direction on the
        orientation sphere, while psi describes the in-plane rotation of the
        particle image. In most biological analyses, preferential orientation
        is primarily associated with rot and tilt, whereas psi often reflects
        rotational variability around the projection axis.

        Defining Angular Regions

        Users define angular intervals that identify the region of orientation
        space targeted for particle removal. The protocol allows specification
        of minimum and maximum values for rot and tilt, thereby selecting a
        region on the Euler sphere corresponding to the preferred views of
        interest.

        In many practical cryo-EM workflows, users first inspect angular
        distribution plots generated during refinement and then define limits
        surrounding the most populated orientation clusters. Narrow limits allow
        focused removal of highly dominant views, whereas broader limits remove
        larger angular regions and produce more aggressive balancing.

        Optional Psi Filtering

        The protocol also provides optional filtering using psi angles. This is
        generally less important for correcting preferential orientation itself,
        because psi mainly represents in-plane image rotation rather than true
        sampling of distinct projection directions. However, psi selection may
        still be useful in specialized analyses where users wish to isolate or
        remove particles exhibiting restricted rotational behavior within the
        image plane.

        In most biological applications, filtering based only on rot and tilt is
        sufficient and biologically more meaningful.

        Particle Removal Strategy

        Once the angular region is defined, the protocol removes a user-defined
        number of particles from within that region. This allows controlled
        reduction of oversampled views without completely eliminating them from
        the dataset.

        The selective nature of the protocol is important biologically because
        completely removing an orientation may introduce new sampling gaps or
        destabilize refinement. Partial reduction often provides a better
        compromise between balancing angular coverage and preserving sufficient
        signal for accurate reconstruction.

        The protocol can also automatically adapt when fewer particles are
        available than requested, ensuring that all eligible particles are
        handled consistently.

        Outputs and Interpretation

        The output consists of a new particle set with the selected particles
        removed. All remaining particles preserve their original metadata and
        alignment information, making the resulting dataset immediately usable
        for downstream refinement, classification, or validation workflows.

        Comparing reconstructions before and after preferential-view removal can
        help users evaluate the impact of angular imbalance on directional
        resolution, map anisotropy, or structural interpretability. In some
        cases, apparent structural features may weaken after balancing the
        dataset, indicating that they were influenced by orientation bias.

        Practical Recommendations

        In routine cryo-EM analysis, it is usually advisable to begin with
        conservative particle removal rather than aggressive filtering. Removing
        a moderate fraction of dominant views often provides sufficient angular
        balancing while preserving reconstruction stability.

        Users should carefully inspect angular distribution plots before
        selecting removal ranges. Overly broad angular limits may excessively
        reduce particle count and degrade reconstruction quality. Likewise,
        removing too many particles from already sparse datasets may reduce
        overall signal-to-noise ratio and negatively affect refinement.

        This protocol is especially valuable for methodological experiments,
        validation studies, and investigations of anisotropic resolution. It is
        less commonly used as a standard preprocessing step in routine
        high-quality datasets with already balanced orientation distributions.

        Final Perspective

        Preferential orientation is fundamentally a biological and biophysical
        phenomenon arising from how particles interact with cryo-EM grids and
        interfaces. The Remove Preferential Views protocol provides a practical
        way to study and partially compensate for these effects by selectively
        reducing dominant orientations. When used carefully, it can help improve
        angular balance, support more isotropic reconstructions, and provide
        deeper insight into the relationship between particle orientation
        distributions and final map quality.
    """
    _label = 'remove preferential views'
    _devStatus = PROD
    _possibleOutputs = outputs

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      pointerCondition='hasAlignmentProj',
                      important=True,
                      label='Input particles',
                      pointerClass='SetOfParticles',
                      help='Provide a set of particles.')
        form.addParam('numToRemove', params.IntParam,
                      default=0,
                      label="Number of particles to remove",
                      help="Number of particles to remove "
                           "WITHIN the limits below.")

        group = form.addGroup("Remove views WITHIN the limits below")
        line = group.addLine('Rot')
        line.addParam('rotMin', params.IntParam, default=-180,
                      label='min')
        line.addParam('rotMax', params.IntParam, default=180,
                      label='max')

        line = group.addLine('Tilt')
        line.addParam('tiltMin', params.IntParam, default=0,
                      label='min')
        line.addParam('tiltMax', params.IntParam, default=180,
                      label='max')

        group.addParam('removePsi', params.BooleanParam, default=False,
                       label="Remove views with specific in-plane rotation?",
                       help="Particle orientation on Euler sphere is "
                            "defined by rot and tilt angles. Psi is for "
                            "in-plane rotation only. Select *Yes* "
                            "if you want to provide psi limits.")
        line = group.addLine('Psi', condition='removePsi')
        line.addParam('psiMin', params.IntParam, default=-180,
                      condition='removePsi', label='min')
        line.addParam('psiMax', params.IntParam, default=180,
                      condition='removePsi', label='max')

    # -------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.processAnglesStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    def convertInputStep(self):
        inputParts = self.inputParticles.get()
        self.rotDict, self.tiltDict, self.psiDict = {}, {}, {}
        for part in inputParts:
            alignment = part.getTransform()
            matrix = np.linalg.inv(alignment.getMatrix())
            angles = -np.rad2deg(tfs.euler_from_matrix(matrix, axes='szyz'))
            self.rotDict[part.getObjId()] = angles[0]
            self.tiltDict[part.getObjId()] = angles[1]
            self.psiDict[part.getObjId()] = angles[2]

    def processAnglesStep(self):
        self.removedList = []
        ptclToRemove = self.numToRemove.get()
        rmin, rmax = self.rotMin.get(), self.rotMax.get()
        tmin, tmax = self.tiltMin.get(), self.tiltMax.get()
        pmin, pmax = self.psiMin.get(), self.psiMax.get()

        for (k, r), (k2, t), (k3, p) in zip(self.rotDict.items(),
                                            self.tiltDict.items(),
                                            self.psiDict.items()):
            if not self.removePsi:
                # check only rot & tilt
                if (self.withinLimits(r, rmin, rmax) and
                        self.withinLimits(t, tmin, tmax)):
                    self.removedList.append(k)
            else:
                if (self.withinLimits(r, rmin, rmax) and
                        self.withinLimits(t, tmin, tmax) and
                        self.withinLimits(p, pmin, pmax)):
                    self.removedList.append(k)

        if ptclToRemove > len(self.removedList):
            self.info("Number to remove (%d) is more than maximum available (%d). "
                      "Removing %d particles..." % (ptclToRemove,
                                                    len(self.removedList),
                                                    len(self.removedList)))
        else:
            self.info("Randomly removing %d particles within "
                      "specified limits..." % ptclToRemove)
            self.removedList = sample(self.removedList, ptclToRemove)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        imgSet.setAlignmentProj()
        outImgSet.copyItems(imgSet, updateItemCallback=self._removeViews)
        self._defineOutputs(**{outputs.outputParticles.name: outImgSet})
        self._defineTransformRelation(self.inputParticles, outImgSet)

    def _removeViews(self, item, row):
        if item.getObjId() in self.removedList:
            setattr(item, "_appendItem", False)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, "outputParticles"):
            summary.append("Input particles: %d" % self.inputParticles.get().getSize())
            summary.append("Output particles: %d" % self.outputParticles.get().getSize())

        return summary

    def _validate(self):
        errors = []

        if self.numToRemove > self.inputParticles.get().getSize():
            errors.append("You cannot remove more particles "
                          "than provided in the input.")

        if (self.rotMin < -180 or self.psiMin < -180 or
                self.rotMax > 180 or self.tiltMax > 180 or
                self.psiMax > 180 or self.tiltMin < 0):
            errors.append("Angles must be within the following limits:\n\n"
                          "-180 < rot < 180\n0 < tilt < 180\n-180 < psi < 180\n")

        return errors

    # -------------------------- UTILS functions ------------------------------
    def withinLimits(self, value, min, max):
        return min < float(value) < max
