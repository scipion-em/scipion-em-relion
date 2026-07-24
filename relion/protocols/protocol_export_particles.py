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

import os

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD
from pwem.constants import ALIGN_NONE
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtProcessParticles

import relion.convert as convert
from ..constants import STACK_MULT, STACK_ONE
from .protocol_base import ProtRelionBase


class ProtRelionExportParticles(ProtProcessParticles, ProtRelionBase):
    """
    Exports particle datasets from Scipion into RELION-compatible STAR
    metadata and optional particle stack files for external processing,
    sharing, or archival purposes.

    AI Generated:

    Export Particles (ProtRelionExportParticles) — User Manual
        Overview

        The Export Particles protocol prepares cryo-EM particle datasets
        for use outside the Scipion environment by converting them into
        RELION-compatible STAR metadata and associated image stacks.
        Its main purpose is to facilitate interoperability between
        software packages, computational facilities, collaborators, and
        long-term storage workflows.

        In practical cryo-EM research, particle export is commonly
        required when transferring datasets between institutions,
        continuing processing in external RELION installations,
        performing specialized downstream analyses, or preparing data
        for publication and public deposition. The protocol ensures
        that particle metadata, imaging parameters, and optional
        alignment information remain organized in a standardized and
        portable format.

        Inputs and General Workflow

        The protocol requires a particle dataset that has already been
        generated within a Scipion workflow. These particles may
        originate from extraction, classification, refinement, or
        polishing stages and can contain additional metadata such as
        alignments, CTF parameters, optics information, and acquisition
        settings.

        During export, the protocol generates STAR metadata files and,
        depending on user preference, may also generate particle stack
        files in MRCS format. The resulting exported dataset can then
        be imported directly into RELION or used by compatible cryo-EM
        software environments.

        Alignment Information

        The protocol optionally preserves alignment information during
        export. This capability is biologically important because
        alignment parameters encode the estimated orientations and
        translations of particles relative to reconstructed structures.

        Including alignment information is typically recommended when
        exporting particles for continued refinement, focused
        classification, heterogeneous analysis, or structural
        interpretation. Preserving these parameters allows downstream
        workflows to continue from previously optimized orientations
        rather than restarting alignment procedures from the beginning.

        In contrast, users may choose to omit alignment information
        when preparing datasets for independent reprocessing,
        benchmarking studies, or unbiased exploratory analysis.

        Particle Stack Organization

        The protocol supports several strategies for organizing binary
        particle stacks. Users may export only metadata, generate a
        single consolidated stack, or create multiple stacks organized
        according to acquisition grouping.

        A single stack simplifies file management and is often
        convenient for transfer between systems or smaller projects.
        This approach is particularly useful when datasets are intended
        for straightforward downstream processing or archival storage.

        Multiple stack organization may be preferable for large-scale
        cryo-EM projects, streaming environments, or facility-based
        workflows where separating particles by micrograph or
        acquisition batch improves scalability and data management.

        Metadata Consistency and Interoperability

        A major advantage of the export workflow is the preservation of
        metadata consistency across software environments. Correct
        handling of pixel size, optics groups, alignment parameters,
        and acquisition settings is essential for maintaining
        reproducibility in cryo-EM processing.

        In modern high-resolution workflows, subtle metadata
        inconsistencies can lead to refinement instability, inaccurate
        scaling, or resolution loss. Proper export therefore plays an
        important role in ensuring that downstream analyses remain
        scientifically reliable.

        Outputs and Their Interpretation

        The protocol produces one or more STAR files describing the
        particle dataset together with optional binary image stacks.
        These exported files contain the information required for
        downstream RELION processing and external data exchange.

        Biologically, the exported dataset represents a portable
        snapshot of the particle processing stage at the moment of
        export. Depending on the workflow, this may correspond to raw
        extracted particles, cleaned datasets after classification, or
        highly refined particles prepared for advanced reconstruction.

        Practical Recommendations

        In most routine cryo-EM workflows, preserving alignment
        information is advisable when the goal is to continue
        refinement or structural analysis outside Scipion. However,
        when preparing datasets for independent validation or
        benchmarking, exporting particles without alignments may reduce
        potential processing bias.

        For small and medium-sized projects, a single particle stack is
        often easier to manage and transfer. In contrast, very large
        datasets may benefit from multiple-stack organization to avoid
        extremely large binary files and improve compatibility with
        distributed computing environments.

        Before sharing exported datasets with collaborators or
        depositing them in repositories, users should verify that pixel
        sizes, particle counts, and metadata consistency remain correct
        after export.

        Collaborative and Facility Workflows

        The protocol is especially useful in collaborative cryo-EM
        environments where datasets move between computational
        facilities, laboratories, and software ecosystems. Exported
        particle packages provide a standardized mechanism for
        exchanging data while preserving essential reconstruction
        metadata.

        In facility pipelines, this export capability also supports
        reproducibility by enabling researchers to archive the precise
        particle dataset associated with a published reconstruction or
        processing milestone.

        Final Perspective

        Particle export is an essential interoperability step in modern
        cryo-EM analysis. Reliable preservation of metadata, alignment
        information, and particle organization ensures that structural
        interpretation and downstream refinement can proceed accurately
        across different computational environments. Careful export
        practices contribute directly to reproducibility, collaboration,
        and long-term scientific reliability.
    """

    _label = 'export particles'
    _devStatus = PROD
    PTCLS_STAR_FILE = 'particles_%06d.star'

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        form.addSection(label='Input')

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Select the input images from the project.')

        form.addParam('useAlignment', params.BooleanParam, default=True,
                      label='Write alignment information?',
                      help='If *Yes* the alignment information (2D or 3D) '
                           'will be written to the resulting .star file if '
                           'the particles contains such information.')

        form.addParam('stackType', params.EnumParam,
                      choices=["Don't write stacks",
                               "Write multiple stacks",
                               "Write a single stack"], default=STACK_MULT,
                      display=params.EnumParam.DISPLAY_LIST,
                      label="Binary stack files",
                      help="If *Don't write stacks* is chosen, only the star "
                           "files will be written out. Alternatively, you can "
                           "select to write images into a single stack file or"
                           " several stacks (one per micrograph). ")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        objId = self.inputParticles.get().getObjId()
        self._insertFunctionStep(self.exportParticlesStep, objId, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def exportParticlesStep(self, particlesId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        pwutils.cleanPath(self._getExportPath())
        pwutils.makePath(self._getExportPath())
        imgSet = self.inputParticles.get()
        self._stackType = self.stackType.get()
        self._ih = ImageHandler()
        self._stackDict = {}

        alignType = imgSet.getAlignment() if self.useAlignment else ALIGN_NONE
        outputDir = None
        outputStack = None
        postprocessImageRow = None

        if self._stackType == STACK_ONE:
            outputStack = self._getExportPath('Particles/particles.mrcs')
            pwutils.makePath(self._getExportPath("Particles"))

        elif self._stackType == STACK_MULT:
            postprocessImageRow = self._postprocessImageRow
            outputDir = self._getExportPath("Particles")

        # Create links to binary files and write the relion .star file
        convert.writeSetOfParticles(
            imgSet, self._getStarFile(),
            outputDir=outputDir,
            outputStack=outputStack,
            alignType=alignType,
            postprocessImageRow=postprocessImageRow,
            fillMagnification=True,
            forceConvert=True)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []
        return validateMsgs

    def _summary(self):
        summary = []

        if os.path.exists(self._getStarFile()):
            summary.append("Output is written to: \n%s\n" %
                           os.path.abspath(self._getExportPath()))
            summary.append("Pixel size: *%0.3f*" % self._getPixelSize())
        else:
            summary.append("No output generated yet.")

        return summary

    # --------------------------- UTILS functions -----------------------------
    def _postprocessImageRow(self, img, row):
        """ Stack fn should be relative to Export.
        Only relevant when saving multiple stacks. """
        convert.relativeFromFileName(row, self._getExportPath())

    def _getExportPath(self, *paths):
        return os.path.join(self._getPath('Export'), *paths)

    def _getStarFile(self):
        return self._getExportPath(self.PTCLS_STAR_FILE % self.getObjId())

    def _getPixelSize(self):
        return self.inputParticles.get().getSamplingRate()
