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

import os

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD

import relion.convert as convert
from .protocol_base import ProtRelionBase


class ProtRelionExportCoordinates(ProtRelionBase):
    """
    Exports particle coordinates into RELION-compatible STAR files so they can
    be reused outside the Scipion environment in downstream cryo-EM workflows.

    AI Generated:

    Export Coordinates (ProtRelionExportCoordinates) — User Manual
        Overview

        The Export Coordinates protocol is designed to transfer particle picking
        information from a Scipion project into a format that can be directly
        understood by RELION and related cryo-EM software ecosystems. Its main
        purpose is to preserve coordinate information associated with
        micrographs so that particle extraction, reprocessing, or independent
        analysis can continue outside the original workflow environment.

        In practical biological workflows, coordinate export is commonly used
        when particle picking has been completed in one software platform but
        subsequent processing will be performed elsewhere. This situation often
        arises in collaborative projects, facility-based processing pipelines,
        or comparative analyses where multiple cryo-EM software packages are
        combined to exploit their complementary strengths.

        Inputs and General Workflow

        The protocol requires a set of particle coordinates associated with one
        or more micrographs. These coordinates usually represent particle
        positions identified either manually or through automated particle
        picking procedures. The export process preserves the spatial location
        of each selected particle relative to its original micrograph.

        During execution, the protocol generates STAR files organized in a way
        that is compatible with RELION conventions. Each micrograph receives an
        associated coordinate file, allowing downstream software to identify
        particle positions consistently and reproducibly. This organization is
        especially important in large cryo-EM datasets where thousands of
        micrographs and millions of particles may be involved.

        Biological Context and Use Cases

        Coordinate export plays an important role in maintaining reproducibility
        across cryo-EM workflows. Particle coordinates define the initial set
        of candidate projections used for all later processing steps, including
        extraction, classification, refinement, and reconstruction. Consistency
        in coordinate handling is therefore essential to ensure that results
        obtained in different software environments remain biologically
        comparable.

        Researchers frequently use exported coordinates when testing different
        extraction strategies, comparing particle-picking algorithms, or
        migrating projects between computational infrastructures. In facility
        settings, coordinate export also facilitates data sharing between users,
        processing teams, and external collaborators.

        Practical Recommendations

        Before exporting coordinates, users should verify that the associated
        micrographs are correctly linked and consistently named. Stable file
        organization greatly simplifies downstream processing and avoids errors
        when importing the exported STAR files into other software packages.

        It is also advisable to confirm that particle picking quality has been
        visually validated before export. Since all subsequent cryo-EM analyses
        depend on the accuracy of the selected coordinates, poor particle
        selection at this stage can negatively affect classification quality,
        reconstruction resolution, and biological interpretation.

        When working with very large datasets, maintaining consistent naming
        conventions between micrographs and coordinate files becomes especially
        important. Proper organization ensures compatibility with automated
        processing pipelines and reduces the likelihood of mismatches during
        particle extraction.

        Outputs and Interpretation

        The protocol produces a collection of RELION-compatible STAR coordinate
        files organized for external use. These outputs can be imported into
        RELION or other compatible cryo-EM software for particle extraction or
        further image processing.

        The exported coordinates preserve the relationship between particles
        and their originating micrographs, ensuring that downstream workflows
        retain the spatial context required for accurate extraction and
        reconstruction.

        Final Perspective

        In cryo-EM image processing, coordinate export represents an important
        interoperability step that allows particle-picking information to move
        efficiently between software environments. Reliable coordinate transfer
        supports reproducibility, collaborative processing, and flexible
        workflow design, all of which are increasingly important in modern
        structural biology projects.
    """

    _label = 'export coordinates'
    _devStatus = PROD

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label="Input coordinates",
                      help='Select the SetOfCoordinates ')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        objId = self.inputCoordinates.get().getObjId()
        self._insertFunctionStep(self.exportCoordsStep, objId, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def exportCoordsStep(self, coordsId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        """
        pwutils.cleanPath(self._getExportPath())
        pwutils.makePath(self._getExportPath())

        convert.writeSetOfCoordinates(self._getExportPath(),
                                      self.getCoords(),
                                      self._getMicPos)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        return validateMsgs

    def _summary(self):
        summary = []

        summary.append("Output is written to: \n%s\n" %
                       os.path.abspath(self._getExportPath()))

        return summary
    
    # --------------------------- UTILS functions -----------------------------
    def _getMicPos(self, mic):
        fileName = pwutils.removeBaseExt(mic.getFileName()) + "_coords.star"
        return fileName

    def getCoords(self):
        return self.inputCoordinates.get()

    def _getExportPath(self, *paths):
        return os.path.join(self._getPath('Export'), *paths)
