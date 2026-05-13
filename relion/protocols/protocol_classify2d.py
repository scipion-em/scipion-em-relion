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

from pyworkflow.constants import PROD
from pwem.constants import ALIGN_2D
from pwem.objects import SetOfClasses2D
from pwem.protocols import ProtClassify2D

import relion.convert as convert
from .protocol_base import ProtRelionBase


class outputs(Enum):
    outputClasses = SetOfClasses2D


class ProtRelionClassify2D(ProtRelionBase, ProtClassify2D):
    """
    Performs unsupervised 2D classification of cryo-EM particle images using
    the Bayesian refinement framework implemented in Relion. The protocol
    groups similar particle projections into representative 2D classes that
    capture common structural views, improve signal-to-noise ratio, and help
    identify heterogeneity, contaminants, or damaged particles within a
    dataset.

    AI Generated:

    Relion 2D Classification (ProtRelionClassify2D) — User Manual
        Overview

        The Relion 2D Classification protocol organizes particle images into
        statistically consistent groups based on similarity of their projected
        structural features. In cryo-EM workflows, this step is one of the
        most important stages for evaluating dataset quality because it allows
        users to separate well-defined particle views from noise, ice
        contamination, aggregation, or incorrectly picked particles.

        The protocol uses the Bayesian optimization strategy developed in
        Relion to estimate both particle alignment parameters and class
        assignments simultaneously. Rather than relying on rigid deterministic
        classification, the method evaluates probabilities for each particle
        belonging to each class, resulting in stable and biologically
        meaningful averages even in noisy experimental datasets.

        Biological Purpose and Typical Applications

        In practical cryo-EM analysis, 2D classification is commonly used as
        an early quality-control stage after particle extraction. Well-defined
        class averages reveal whether the particles contain recognizable
        structural information and whether the dataset is suitable for further
        refinement or reconstruction.

        Biological users frequently employ this protocol to remove false
        positives from automatic picking, eliminate damaged particles, and
        isolate subsets corresponding to different orientations of the same
        macromolecule. In many projects, this step also provides the first
        visual confirmation of structural integrity before computationally
        intensive 3D analysis begins.

        For highly heterogeneous samples, 2D classification may additionally
        reveal distinct particle populations, preferred orientations, or
        flexibility. Although the protocol does not directly reconstruct
        three-dimensional maps, the quality and diversity of the resulting
        classes strongly influence all downstream stages of cryo-EM processing.

        Inputs and Dataset Considerations

        The protocol requires a set of extracted particle images as input.
        Optimal results are usually obtained when particles have already been
        motion-corrected, contrast-transfer-function estimated, and reasonably
        centered during extraction. Excessively noisy particles or particles
        with very large positional errors may reduce classification quality.

        The number of requested classes should reflect the expected complexity
        and size of the dataset. Smaller numbers of classes generally produce
        smoother and more stable averages, while larger numbers may reveal
        subtle structural variability or rare particle orientations. However,
        requesting too many classes can fragment the dataset and produce noisy
        or poorly populated averages.

        Alignment and Angular Sampling

        During classification, the protocol can optimize translational and
        rotational alignment parameters for each particle. This alignment step
        is biologically important because particles in experimental cryo-EM
        images appear with arbitrary in-plane rotations and small positional
        shifts.

        The translational search range determines how far particles may move
        during alignment. Wider ranges are useful when extracted particles are
        poorly centered, although they increase computational cost. Narrower
        ranges are generally preferred once particle coordinates are already
        accurate.

        Angular sampling controls the precision of in-plane rotational
        searches. Coarse angular sampling accelerates processing and is often
        sufficient during exploratory classification, whereas finer angular
        sampling improves alignment precision and class sharpness at the cost
        of additional computation.

        In some workflows, users may disable image alignment entirely. This is
        primarily useful for specialized analyses where particles are already
        aligned externally or when the goal is to evaluate structural
        variability without introducing additional alignment refinement.

        Gradient-Based Optimization

        The protocol may optionally employ gradient-based optimization methods
        for accelerated refinement. These approaches can substantially reduce
        computation time for large datasets while maintaining classification
        quality. They are especially valuable in facility-scale workflows or
        during rapid exploratory processing.

        However, accelerated optimization strategies should still be evaluated
        critically by biological users. In difficult datasets with strong
        heterogeneity, low signal, or severe preferred orientation,
        conservative refinement strategies may provide more stable and
        interpretable class averages.

        Iterative Refinement and Continuation

        Classification proceeds iteratively, progressively refining particle
        alignments and class averages over multiple cycles. Early iterations
        typically produce broad low-resolution groupings, while later
        iterations sharpen structural details and improve separation between
        classes.

        The protocol supports continuation from previous runs, allowing users
        to extend refinement, increase iterations, or further optimize
        challenging datasets without restarting the entire workflow. This is
        particularly useful for large cryo-EM projects where classification
        may evolve alongside improvements in particle cleaning or extraction.

        Outputs and Biological Interpretation

        The primary output consists of a set of 2D class averages together
        with the particles assigned to each class. These classes summarize the
        dominant structural views present in the dataset and frequently serve
        as the basis for particle selection before 3D reconstruction.

        Well-resolved classes generally display recognizable secondary
        structure features, consistent particle boundaries, and low background
        noise. Poorly resolved or noisy classes often correspond to damaged
        particles, contaminants, overlapping particles, or alignment failure.

        Biological interpretation should consider both the visual quality and
        the distribution of particles across classes. Large populations of
        high-quality classes usually indicate a healthy dataset, whereas
        excessive fragmentation or dominance of junk classes may suggest
        problems in sample preparation, imaging conditions, or particle
        picking.

        Practical Recommendations

        For routine cryo-EM processing, it is generally advisable to begin
        with a moderate number of classes and standard alignment settings.
        After inspecting the resulting averages, users may refine the analysis
        by increasing class number, tightening alignment parameters, or
        removing poor particles before additional rounds of classification.

        Datasets with severe heterogeneity often benefit from multiple rounds
        of classification, where obvious contaminants are removed first and
        subsequent runs focus on improving biologically relevant subsets.
        Careful visual inspection remains essential throughout the process
        because automated classification alone cannot fully determine the
        biological relevance of each class.

        Final Perspective

        For most cryo-EM workflows, 2D classification represents a critical
        bridge between raw particle extraction and high-resolution structural
        analysis. Beyond its computational role, it provides direct biological
        insight into particle quality, conformational diversity, and dataset
        integrity. Thoughtful interpretation of class averages and careful
        selection of high-quality particles are essential for achieving robust
        and biologically meaningful downstream reconstructions.
    """

    _label = '2D classification'
    _devStatus = PROD
    _possibleOutputs = outputs
    IS_2D = True
    OUTPUT_TYPE = SetOfClasses2D
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
        self.ClassFnTemplate = '%(ref)03d@%(rootDir)s/relion_it%(iter)03d_classes.mrcs'

    # --------------------------- INSERT steps functions ----------------------
    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        if self.doImageAlignment:
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step'] = self.offsetSearchStepPix.get() * self._getSamplingFactor()
            args['--psi_step'] = self.inplaneAngularSamplingDeg.get() * self._getSamplingFactor()

            if self.allowCoarserSampling:
                args['--allow_coarser_sampling'] = ''

        else:
            args['--skip_align'] = ''

    # --------------------------- STEPS functions -----------------------------
    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses2D from a given iteration. """
        classLoader = convert.ClassesLoader(self, ALIGN_2D)
        classLoader.fillClassesFromIter(clsSet, iteration)

    def createOutputStep(self):
        classes2D = self._createSetOfClasses2D(self.inputParticles)
        self._fillClassesFromIter(classes2D, self._lastIter())
        
        self._defineOutputs(**{outputs.outputClasses.name: classes2D})
        self._defineSourceRelation(self.inputParticles, classes2D)
        
    # --------------------------- INFO functions ------------------------------
    def _validateNormal(self):
        errors = []
        if self.useGradientAlg and self.numberOfMpi > 1:
            errors.append("Gradient refinement (running the VDAM algorithm) "
                          "is not supported together with MPI")

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
        summary = [
            "Input Particles: %s" % self.getObjectTag('inputParticles'),
            "Classified into *%d* classes." % self.numberOfClasses,
            "Output set: %s" % self.getObjectTag('outputClasses')
            ]
        
        return summary
    
    def _summaryContinue(self):
        summary = ["Continue from iteration %01d" % self._getContinueIter()]
        return summary
    
    def _methods(self):
        methods = ''
        if hasattr(self, 'outputClasses'):
            methods += "We classified input particles %s (%d items) " % (
                self.getObjectTag('inputParticles'),
                self.inputParticles.get().getSize())
            methods += "into %d classes using Relion Classify2d. " % self.numberOfClasses
            methods += 'Output classes: %s' % self.getObjectTag('outputClasses')
        return [methods]
