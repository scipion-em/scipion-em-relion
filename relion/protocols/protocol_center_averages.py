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
    Aligns two-dimensional class averages by centering them according to
    their center of mass. This procedure is commonly used in cryo-EM image
    processing workflows to improve the visual consistency and positional
    normalization of class averages before downstream analysis, refinement,
    visualization, or comparison. By repositioning the particle signal to a
    common central location, the protocol facilitates more reliable
    interpretation of structural features and helps reduce variability caused
    by translational offsets.

    AI Generated:

    Center Averages (ProtRelionCenterAverages) - User Manual
        Overview

        The Center Averages protocol recenters a set of two-dimensional class
        averages using the center-of-mass alignment strategy implemented in
        Relion image processing tools. The main objective is to place the
        dominant particle density at the center of each image so that all
        averages share a more consistent spatial reference frame.

        In practical cryo-EM workflows, class averages are often generated
        after several rounds of particle alignment and classification.
        Although these averages usually represent coherent structural views,
        they may still contain small translational displacements caused by
        image variability, alignment uncertainty, or heterogeneity in the
        dataset. Centering the averages improves visual uniformity and makes
        subsequent interpretation more straightforward.

        Biological Motivation

        From a biological perspective, centered class averages are easier to
        inspect and compare. Structural motifs, domains, and conformational
        features become more visually consistent across the dataset when the
        particle occupies a stable position within the image frame. This is
        particularly useful during exploratory analysis, quality assessment,
        and preparation of publication-quality figures.

        In many workflows, centered averages are also beneficial before
        downstream computational steps such as template generation, initial
        model estimation, classification refinement, or particle selection.
        Ensuring that averages share a common positional reference can improve
        the robustness of later alignment stages.

        Inputs and Expected Data

        The protocol requires a set of two-dimensional class averages as
        input. These averages are typically generated from previous
        classification procedures and should already represent meaningful
        particle views. The protocol is intended for averaged particle images
        rather than raw micrographs or individual particles.

        The quality of the input averages strongly influences the usefulness
        of the centering operation. Well-defined averages with compact signal
        regions are generally centered reliably, whereas highly noisy or
        heterogeneous averages may produce less meaningful repositioning.

        General Workflow

        During execution, each average is analyzed to determine the spatial
        distribution of its intensity. The protocol estimates the center of
        mass of the particle signal and shifts the image so that the detected
        center is relocated toward the middle of the image frame.

        The procedure is designed to preserve the identity and appearance of
        each average while improving positional consistency across the entire
        dataset. No classification, averaging, or structural modification is
        introduced during this process.

        Outputs and Interpretation

        The protocol produces a new set of centered class averages that retain
        the same structural information as the input set but with improved
        spatial alignment. The resulting averages can be used directly for
        visualization, further refinement, template creation, or downstream
        cryo-EM processing tasks.

        Biologically, the centered averages should not be interpreted as new
        structural states or altered reconstructions. The operation only
        standardizes image positioning and does not modify the underlying
        particle information.

        Practical Considerations

        In most biological applications, centering class averages is a simple
        but valuable preprocessing step that improves dataset consistency and
        interpretability. It is particularly useful when averages appear
        displaced within the image box or when multiple averages need to be
        visually compared side by side.

        Users should nevertheless inspect the outputs visually, especially in
        datasets containing elongated, asymmetric, flexible, or fragmented
        particles. In such situations, the estimated center of mass may not
        always correspond to the biologically most informative alignment
        position.

        Final Perspective

        For cryo-EM practitioners, centering class averages is often a small
        but important refinement step that improves the clarity and coherence
        of downstream analyses. By placing particle projections into a common
        positional frame, the protocol contributes to cleaner visualization,
        more reliable comparison between classes, and smoother integration
        into broader single-particle analysis workflows.
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
                                 self.inputAverages.get().getObjId(),
                                 needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

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
