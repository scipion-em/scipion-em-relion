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
from emtable import Table

from pyworkflow.constants import PROD
from pwem.constants import ALIGN_PROJ
from pwem.objects import SetOfClasses3D, SetOfVolumes
from pwem.protocols import ProtClassify3D

import relion.convert as convert
from .protocol_base import ProtRelionBase


class outputs(Enum):
    outputClasses = SetOfClasses3D
    outputVolumes = SetOfVolumes


class ProtRelionClassify3D(ProtClassify3D, ProtRelionBase):
    """
    Performs 3D classification of cryo-EM particle images using the
    Bayesian refinement framework implemented in RELION. The protocol
    separates heterogeneous particle populations into distinct 3D
    structural classes, allowing the identification of different
    conformational, compositional, or quality-related states within a
    dataset.

    AI Generated:

    3D Classification (ProtRelionClassify3D) — User Manual
        Overview

        The 3D Classification protocol applies RELION Bayesian methods
        to classify cryo-EM particles into multiple three-dimensional
        structural groups. Its primary purpose is to identify structural
        variability within experimental datasets and to separate
        particles that represent different conformations, assemblies,
        or data quality levels.

        In cryo-EM workflows, biological samples frequently contain a
        mixture of structural states rather than a single homogeneous
        population. Different ligand occupancies, flexible domains,
        assembly intermediates, or damaged particles may coexist in the
        same dataset. This protocol helps disentangle those mixed
        populations by assigning particles to distinct 3D classes based
        on their statistical agreement with evolving reference volumes.

        For many biological projects, 3D classification represents one
        of the most important decision-making stages in the entire
        workflow because it determines which subsets of particles are
        suitable for high-resolution refinement and which subsets
        correspond to alternative or low-quality states.

        Inputs and General Workflow

        The protocol requires a set of input particles together with an
        initial reference volume. Using Bayesian optimization, the
        particles are iteratively aligned and classified into multiple
        3D classes that progressively improve during refinement.

        The resulting classes can represent biologically meaningful
        conformations, distinct assembly compositions, flexible
        structural arrangements, or experimental artifacts. Because the
        classification process is probabilistic, particles are assigned
        according to their likelihood of belonging to each structural
        group rather than through rigid deterministic rules.

        In practical cryo-EM workflows, the quality of the initial
        reference strongly influences convergence behavior. References
        should ideally represent the overall architecture of the
        particle while avoiding strong model bias. Excessively detailed
        references may artificially favor particular conformations and
        reduce the ability to discover unexpected structural states.

        Image Alignment and Orientation Search

        The protocol can perform simultaneous image alignment during
        classification, allowing particle orientations and translations
        to be refined together with class assignments. This is the most
        common strategy in biological workflows because particles often
        require continuous orientation optimization while heterogeneity
        is being resolved.

        Angular and translational searches determine how broadly the
        protocol explores possible orientations and shifts. Wider
        searches improve robustness when particle orientations are
        uncertain, although they increase computational cost. Narrower
        searches are faster and more precise when approximate alignment
        is already available.

        Local angular searches may be particularly useful in later
        stages of classification when particles are already reasonably
        aligned. In such cases, restricting orientation exploration can
        stabilize refinement and improve convergence toward subtle
        conformational differences.

        Biological Interpretation of Classes

        The interpretation of resulting classes requires careful
        biological judgment. Some classes may represent true structural
        states, while others may correspond to damaged particles,
        preferred orientations, contamination, or partially aligned
        subsets.

        In many cryo-EM projects, the protocol is used iteratively.
        Initial classifications often separate broad structural groups,
        while subsequent focused classifications refine specific
        conformational differences or isolate minor populations of
        biological interest.

        Flexible complexes, membrane proteins, and dynamic assemblies
        particularly benefit from 3D classification because averaging
        all particles together may obscure meaningful structural
        variability. By separating heterogeneous states, the protocol
        enables more accurate structural interpretation and higher
        quality downstream refinement.

        Outputs and Their Interpretation

        The protocol produces a set of classified particles together
        with representative 3D volumes corresponding to each structural
        class. Each volume reflects the average structural information
        of the particles assigned to that class.

        The resulting classes can be inspected visually to identify
        biologically meaningful conformations, assess reconstruction
        quality, and determine which particle subsets should continue
        toward high-resolution refinement.

        Resolution estimates associated with the classifications help
        evaluate convergence and structural consistency. However,
        biological interpretation should not rely solely on nominal
        resolution values. Visual inspection of structural features and
        consistency across independent refinements remain essential.

        Continuation and Iterative Refinement

        The protocol supports continuation from previous classification
        iterations, allowing users to extend or refine existing
        classification runs. This capability is valuable when exploring
        difficult heterogeneous datasets that require gradual
        optimization.

        Iterative classification strategies are common in cryo-EM
        analysis. Researchers often begin with broad classifications to
        remove damaged or irrelevant particles and later perform more
        focused classifications on selected subsets to resolve subtle
        conformational variability.

        Practical Recommendations

        In routine biological practice, selecting an appropriate number
        of classes is one of the most important decisions. Too few
        classes may merge distinct conformations, whereas too many
        classes can fragment particle populations and reduce the signal
        available for reconstruction.

        For highly heterogeneous datasets, beginning with a relatively
        large number of classes is often useful for exploratory
        analysis. Classes can later be merged or refined depending on
        their biological relevance and reconstruction quality.

        Care should also be taken when interpreting poorly populated
        classes. Small classes may reveal rare but meaningful
        conformations, but they may also arise from noise, alignment
        instability, or imaging artifacts.

        Final Perspective

        Three-dimensional classification is one of the central tools of
        modern cryo-EM because it transforms heterogeneous particle
        datasets into interpretable structural populations. Beyond its
        computational role, it provides direct biological insight into
        molecular flexibility, compositional variability, and dynamic
        structural behavior. Careful interpretation of the resulting
        classes, combined with iterative refinement and visual
        validation, is essential for obtaining biologically reliable
        conclusions from cryo-EM experiments.
    """

    _label = '3D classification'
    _devStatus = PROD
    _possibleOutputs = outputs
    CHANGE_LABELS = ['rlnChangesOptimalOrientations',
                     'rlnChangesOptimalOffsets',
                     'rlnOverallAccuracyRotations',
                     'rlnOverallAccuracyTranslationsAngst',
                     'rlnChangesOptimalClasses']
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
    
    # -------------------------- INSERT steps functions -----------------------
    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        if self.doImageAlignment:
            args['--healpix_order'] = self.angularSamplingDeg.get()
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step'] = self.offsetSearchStepPix.get() * self._getSamplingFactor()

            # check if sigma_ang is in extra params
            # before adding the default value
            if self.localAngularSearch:
                if self.relaxSymm.get():
                    args['--relax_sym'] = self.relaxSymm.get()
                if self.extraParams.hasValue():
                    if self.extraParams.get().find("--sigma_ang") == -1:
                        args['--sigma_ang'] = self.localAngularSearchRange.get() / 3.
                else:
                    args['--sigma_ang'] = self.localAngularSearchRange.get() / 3.

            if self.allowCoarserSampling:
                args['--allow_coarser_sampling'] = ''

        else:
            args['--skip_align'] = ''
    
    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        partSet = self.inputParticles
        classes3D = self._createSetOfClasses3D(partSet)
        self._fillClassesFromIter(classes3D, self._lastIter())
        
        self._defineOutputs(**{outputs.outputClasses.name: classes3D})
        self._defineSourceRelation(partSet, classes3D)

        # create a SetOfVolumes and define its relations
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(partSet.get().getSamplingRate())
        
        for class3D in classes3D:
            vol = class3D.getRepresentative()
            vol.setObjId(class3D.getObjId())
            volumes.append(vol)
        
        self._defineOutputs(**{outputs.outputVolumes.name: volumes})
        self._defineSourceRelation(partSet, volumes)
        
        if not self.doContinue:
            self._defineSourceRelation(self.referenceVolume, classes3D)
            self._defineSourceRelation(self.referenceVolume, volumes)
    
    # -------------------------- INFO functions -------------------------------
    def _validateNormal(self):
        errors = []
        return errors
    
    def _validateContinue(self):
        errors = []
        continueRun = self.continueRun.get()
        continueRun._initialize()
        lastIter = continueRun._lastIter()
        
        if self.continueIter.get() == 'last':
            continueIter = lastIter
        else:
            continueIter = int(self.continueIter.get())
        
        if continueIter > lastIter:
            errors += ["You can continue only from the iteration %01d or less" % lastIter]
        
        return errors
    
    def _summaryNormal(self):
        summary = []
        it = self._lastIter() or -1
        if it >= 1:
            table = Table(fileName=self._getFileName('model', iter=it),
                          tableName='model_general')
            row = table[0]
            resol = float(row.rlnCurrentResolution)
            summary.append("Current resolution: *%0.2f A*" % resol)

        inputParts = self.inputParticles.get()
        sizeStr = 'None' if inputParts is None else inputParts.getSize()
        summary.append("Input Particles: *%s*\n"
                       "Classified into *%d* 3D classes\n"
                       % (sizeStr, self.numberOfClasses))
        
        return summary
    
    def _summaryContinue(self):
        summary = list()
        summary.append("Continue from iteration %01d" % self._getContinueIter())
        return summary
    
    def _methods(self):
        strline = ''
        if hasattr(self, 'outputClasses'):
            strline += 'We classified %d particles into %d 3D classes using Relion Classify3d. ' %\
                           (self.inputParticles.get().getSize(), self.numberOfClasses.get())
        return [strline]
    
    # -------------------------- UTILS functions ------------------------------
    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses3D from a given iteration. """
        classLoader = convert.ClassesLoader(self, ALIGN_PROJ)
        classLoader.fillClassesFromIter(clsSet, iteration)
