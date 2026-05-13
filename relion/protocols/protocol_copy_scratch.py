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

import os
import numpy as np
from random import sample
from enum import Enum

import pyworkflow.protocol.params as params
from pyworkflow.constants import BETA
from pwem.protocols import EMProtocol
from pwem.objects import SetOfParticles
import pwem.convert.transformations as tfs

from emtools.utils import Process


class outputs(Enum):
    outputParticles = SetOfParticles


class ProtParticlesToScratch(EMProtocol):
    """
    Transfers particle datasets to a scratch storage location in order
    to improve data accessibility and optimize performance during
    cryo-EM processing workflows. The protocol creates a relocated copy
    of the particle stack while preserving the associated metadata and
    structural relationships required for downstream analysis.

    AI Generated:

    Particles To Scratch (ProtParticlesToScratch) — User Manual
        Overview

        The Particles To Scratch protocol is designed to duplicate
        particle datasets into a dedicated scratch or temporary storage
        area commonly used in high-performance computing environments.
        Its primary purpose is to improve data throughput and reduce
        bottlenecks caused by slow network storage during demanding
        cryo-EM processing tasks.

        In many cryo-EM facilities and cluster environments, particle
        stacks may initially reside on long-term storage systems that
        prioritize reliability over speed. Although suitable for data
        preservation, these systems can become limiting during intensive
        computational procedures such as refinement, classification, or
        particle polishing. Relocating particle data to scratch storage
        allows faster access and more efficient parallel execution.

        Inputs and General Workflow

        The protocol requires a set of particles as input. These
        particles may originate from extraction, classification,
        refinement, or any intermediate cryo-EM processing stage. The
        protocol creates a corresponding dataset in scratch storage while
        maintaining the original structural organization and particle
        metadata.

        From a biological perspective, the transferred particles remain
        equivalent to the original dataset. The operation does not alter
        particle alignment, image content, contrast properties, or
        metadata interpretation. Instead, the protocol focuses on data
        relocation to improve computational efficiency and workflow
        scalability.

        Scratch Storage in Cryo-EM Workflows

        Scratch storage plays an important role in large-scale cryo-EM
        processing because modern datasets often contain millions of
        particles and terabytes of associated image data. Operations such
        as iterative refinement or large classification jobs repeatedly
        access particle images, making storage performance a critical
        factor in overall execution speed.

        By moving particle data closer to computational resources, the
        protocol reduces file access latency and minimizes contention on
        shared storage systems. This can significantly improve runtime in
        cluster or GPU-based processing environments where many parallel
        tasks access the same particle stack simultaneously.

        Biological users typically encounter the greatest benefits when
        processing very large datasets, performing multiple refinements,
        or working within distributed computing infrastructures.

        Metadata Preservation

        An important aspect of the protocol is the preservation of
        particle metadata and dataset consistency. The relocated particle
        set retains the same experimental and processing information as
        the original input dataset, ensuring compatibility with
        downstream cryo-EM procedures.

        Alignment parameters, acquisition information, particle identity,
        and dataset relationships remain unchanged. As a result, the
        relocated particles can be used interchangeably with the original
        dataset in subsequent processing steps.

        This preservation is especially important in iterative cryo-EM
        workflows where particles may already contain refined alignment
        parameters or classification assignments that must remain
        consistent across processing stages.

        Outputs and Their Interpretation

        The protocol produces a new particle dataset linked to files
        stored in the scratch location. Biologically and computationally,
        the output represents the same particle population as the input,
        but optimized for rapid access during subsequent calculations.

        The output dataset can be directly used in refinement,
        classification, polishing, reconstruction, or any downstream
        cryo-EM procedure requiring particle images.

        Practical Recommendations

        In routine cryo-EM practice, scratch relocation is most useful
        before computationally intensive stages such as high-resolution
        refinement, Bayesian polishing, or large-scale heterogeneous
        classification. Using scratch storage at these stages can improve
        throughput and reduce waiting times in shared computing
        environments.

        Users should ensure that sufficient scratch storage space is
        available before transferring large datasets. Because scratch
        areas are often temporary, they should not replace long-term data
        preservation systems. Important particle datasets and processing
        results should still be archived in permanent storage.

        It is also advisable to verify institutional policies regarding
        scratch cleanup schedules, as some systems automatically remove
        temporary data after a defined period of inactivity.

        Final Perspective

        Efficient data management is a critical component of modern
        cryo-EM processing pipelines. Although transferring particles to
        scratch storage does not alter biological information directly,
        it can substantially improve the efficiency, scalability, and
        responsiveness of computational workflows. Proper use of scratch
        storage therefore contributes to faster and more reliable
        cryo-EM data analysis in large-scale research environments.
    """
    _label = 'particles to scratch'
    _devStatus = BETA
    _possibleOutputs = outputs

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      important=True,
                      label='Input particles',
                      pointerClass='SetOfParticles',
                      help='Provide a set of particles.')

    # -------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    def createOutputStep(self):
        projName = self.getProject().shortName
        logger = Process.Logger()
        imgSet = self.inputParticles.get()
        fnMap = {}
        dirMap = {}

        def _updatePath(item, row):
            fn = item.getFileName()
            dn, base = os.path.split(fn)
            dstDir = os.path.join(os.environ['SCIPION_SCRATCH'], projName, dn)

            if fn not in fnMap:
                fnMap[fn] = True
                if dn not in dirMap:
                    logger.mkdir(dstDir)
                    dirMap[dn] = dstDir
                logger.cp(fn, dstDir)

            item.setFileName(os.path.join(dstDir, base))

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        imgSet.setAlignment(imgSet.getAlignment())
        outImgSet.copyItems(imgSet, updateItemCallback=_updatePath)
        self._defineOutputs(**{outputs.outputParticles.name: outImgSet})
        self._defineTransformRelation(self.inputParticles, outImgSet)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, "outputParticles"):
            summary.append("Input particles: %d" % self.inputParticles.get().getSize())
            summary.append("Output particles: %d" % self.outputParticles.get().getSize())

        return summary

    def _validate(self):
        errors = []
        return errors
