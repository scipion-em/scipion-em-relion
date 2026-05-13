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
    Relion protocol to automatically identify and select high-quality 2D
    class averages from a previous 2D classification workflow.

    AI Generated:

    2D Class Ranker (ProtRelionSelectClasses2D) - User Manual
        Overview

        The 2D Class Ranker protocol is designed to assist cryo-EM users in
        identifying the most reliable and biologically meaningful 2D class
        averages after a RELION 2D classification job. In large datasets,
        manual inspection of hundreds or thousands of classes can become
        time-consuming and subjective. This protocol provides an automated
        strategy to rank classes according to their predicted quality and
        estimated resolution, helping users rapidly separate informative
        particle populations from noise, contaminants, or poorly aligned
        classes.

        In practical cryo-EM workflows, this step is especially valuable
        during early dataset cleaning and quality assessment. By selecting
        only the best classes, users can enrich particle subsets that are
        more suitable for downstream refinement, ab initio reconstruction,
        heterogeneous analysis, or high-resolution structure determination.
        The protocol is intended both for exploratory processing and for
        large-scale automated pipelines where rapid and reproducible class
        evaluation is essential.

        Inputs and General Workflow

        The protocol requires a completed RELION 2D classification run as
        input. From this classification, the protocol evaluates the produced
        class averages and estimates which classes are most likely to
        represent meaningful structural signal. The ranking process combines
        quality prediction with estimated resolution information to guide
        automatic selection.

        The user defines a minimum score threshold that controls how strict
        the automatic selection should be. Higher thresholds generally retain
        only the cleanest and most reproducible classes, while lower
        thresholds include a broader range of particle populations. In
        practice, moderate thresholds are often useful during exploratory
        processing, whereas stricter thresholds are preferred before
        high-resolution refinement.

        The protocol also allows filtering based on estimated resolution.
        This is biologically important because classes with poor resolution
        often correspond to damaged particles, compositional heterogeneity,
        severe flexibility, contamination, or alignment instability.
        Restricting selection to better-resolved classes usually improves
        downstream reconstruction quality.

        Automatic Selection Strategy

        The ranking system attempts to identify classes that contain
        reproducible structural features rather than random noise or
        experimental artifacts. In cryo-EM workflows, good 2D classes often
        display recognizable particle views, clear secondary-structure
        features, and consistent alignment. Poor classes, by contrast, may
        appear blurry, fragmented, aggregated, or dominated by background
        signal.

        Automatic selection is particularly useful for large datasets where
        manual classification becomes impractical. However, biological
        interpretation remains important. Certain biologically relevant
        conformations or rare orientations may receive lower scores simply
        because they are underrepresented. Users should therefore interpret
        automated selection as an efficient filtering strategy rather than a
        complete replacement for expert inspection.

        Thresholds and Selection Criteria

        The minimum threshold parameter controls the balance between purity
        and completeness. Higher thresholds increase confidence in the
        retained classes but may exclude rare or flexible conformations.
        Lower thresholds preserve more diversity but may also retain noisy
        particles.

        The protocol additionally supports retaining a minimum number of
        particles or classes, even when their scores fall below the selected
        threshold. This is useful when users need to preserve dataset size
        for downstream processing or when the dataset contains substantial
        structural variability. In heterogeneous biological systems,
        aggressive filtering can unintentionally eliminate important
        conformational states.

        Selecting a minimum number of particles prioritizes dataset size,
        which may benefit downstream refinement stability. Selecting a
        minimum number of classes instead prioritizes angular diversity and
        conformational coverage. These two strategies should be chosen
        according to the biological objective of the experiment.

        Resolution Considerations

        Estimated class resolution provides an additional layer of quality
        control. Classes with better estimated resolution usually correspond
        to particles that align consistently and share stable structural
        features. Poor-resolution classes often indicate heterogeneity,
        preferred orientation problems, flexibility, or experimental noise.

        From a biological perspective, users should avoid interpreting
        resolution values as absolute measures of structural correctness.
        Certain flexible complexes or membrane proteins may naturally produce
        lower-resolution classes despite containing biologically meaningful
        information. In such cases, a less restrictive resolution cutoff may
        preserve important functional states.

        Outputs and Their Interpretation

        After execution, the protocol produces a selected subset of 2D
        classes together with the associated particles belonging to those
        classes. The resulting particle subset is intended to represent a
        cleaner and more homogeneous population suitable for downstream
        cryo-EM analysis.

        The selected classes retain their associated quality metrics,
        allowing users to inspect predicted scores and estimated resolutions.
        This information can help guide additional manual curation or support
        decisions about subsequent processing strategies.

        The output particles are typically used for ab initio
        reconstruction, 3D classification, refinement, or heterogeneous
        analysis. Cleaner particle subsets often improve convergence,
        increase reconstruction stability, and reduce the influence of
        contaminants or damaged particles.

        Practical Recommendations

        In routine cryo-EM workflows, it is often advisable to begin with a
        moderate threshold and visually inspect the resulting classes before
        applying more aggressive filtering. Automated ranking is highly
        effective for removing obvious noise classes, but manual inspection
        remains valuable for detecting rare conformations or subtle
        structural features.

        For highly heterogeneous samples, users should be cautious with very
        strict thresholds because biologically meaningful variability may be
        discarded. In contrast, for highly homogeneous particles intended
        for high-resolution refinement, stronger filtering often improves
        final map quality.

        When processing very large datasets, the protocol can substantially
        accelerate workflow efficiency by reducing the amount of manual class
        inspection required. Combining automated ranking with expert
        biological interpretation generally provides the best results.

        Final Perspective

        Automatic 2D class selection is an important quality-control step in
        modern cryo-EM processing pipelines. By identifying particle classes
        with the strongest structural signal, the protocol helps users focus
        downstream analysis on the most reliable and informative particle
        populations. Careful adjustment of score thresholds, resolution
        limits, and dataset retention criteria allows the protocol to be
        adapted to a wide range of biological samples and experimental goals.
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
        form.addParam('minResolution', params.FloatParam, default=25,
                      label='Min. resolution',
                      help='Only classes with a estimated resolution '
                           'better than this value will be selected.')
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
        self._insertFunctionStep(self.runSelectStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

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
                row = clsDict.get(locationToRelion(idx, fn), None)
                if row and row.rlnEstimatedResolution < self.minResolution:
                    return row

                return None

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
