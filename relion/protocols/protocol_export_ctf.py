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

from pwem.objects import SetOfMicrographs
from pwem.protocols import EMProtocol
from pyworkflow.constants import PROD
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params

import relion.convert as convert


class ProtRelionExportCtf(EMProtocol):
    """
    Exports contrast transfer function information together with associated
    micrographs into RELION-compatible STAR files for use outside the current
    processing environment. The protocol allows researchers to transfer CTF
    estimations, metadata, and optional power spectrum density references into
    external RELION workflows while preserving the relationship between imaging
    parameters and their corresponding micrographs.

    AI Generated:

    Export CTF (ProtRelionExportCtf) — User Manual
        Overview

        The Export CTF protocol prepares a set of CTF estimations for use in
        external RELION processing pipelines. Its main purpose is to generate
        standardized STAR metadata files that describe the optical properties
        of cryo-EM micrographs together with the corresponding image references.
        This enables interoperability between processing environments and allows
        datasets generated within one workflow to continue processing in RELION
        without loss of essential microscope information.

        In practical cryo-EM projects, CTF export is commonly required when
        users wish to perform downstream classification, refinement, or particle
        extraction in RELION while maintaining previously estimated defocus and
        optical correction parameters. The protocol is particularly useful in
        collaborative facilities, multi-software workflows, and archival data
        preparation.

        Inputs and General Workflow

        The protocol requires a set of CTF estimations as its primary input.
        These estimations are associated with micrographs that define the imaging
        conditions and acquisition geometry. By default, the export uses the
        same micrographs originally employed during CTF estimation, ensuring
        consistency between metadata and image content.

        Alternatively, users may choose to export the CTF information using a
        different set of micrographs. This is especially useful in workflows
        where dose-weighted or otherwise processed micrographs are preferred for
        downstream analysis while still relying on previously computed CTF
        parameters. In these situations, the user should ensure that the
        replacement micrographs correspond correctly to the original acquisition
        geometry and sampling conditions.

        Handling of Power Spectrum Density Information

        The protocol can also preserve associated power spectrum density
        references when they are available. These files provide visual and
        analytical support for evaluating the quality of the CTF estimation and
        are often useful for validation, troubleshooting, or record keeping.

        From a biological and experimental perspective, maintaining access to
        these references can be important when comparing imaging quality across
        datasets, microscope sessions, or acquisition strategies. However,
        absence of these references does not prevent the export of the core CTF
        metadata itself.

        Metadata Consistency and Reliability

        Successful downstream processing depends strongly on maintaining
        consistency between exported micrographs and their associated CTF
        parameters. If micrographs are missing, inaccessible, or incorrectly
        matched, the resulting dataset may become unusable or scientifically
        misleading.

        Biological users should therefore verify that all exported micrographs
        correspond to the same acquisition conditions under which the CTF values
        were estimated. Differences in pixel size, preprocessing strategy, or
        image orientation can negatively affect downstream refinement and
        reconstruction quality.

        Outputs and Their Interpretation

        After execution, the protocol produces RELION-compatible STAR files
        together with references to the exported micrographs and optional power
        spectrum density data. These outputs provide the optical metadata
        required by RELION for subsequent particle processing and reconstruction
        steps.

        The exported STAR files preserve the relationship between each
        micrograph and its associated CTF parameters, allowing external software
        to correctly interpret imaging conditions such as defocus and sampling
        information. The resulting dataset is intended to function as a portable
        representation of the experimental optical metadata.

        Practical Recommendations

        In routine cryo-EM practice, it is generally advisable to export CTF
        information together with the same micrographs used during estimation,
        as this minimizes the risk of metadata inconsistencies. When alternative
        micrographs are required, users should carefully confirm that image
        dimensions, pixel sizes, and acquisition identifiers remain compatible.

        Retaining power spectrum density references is recommended whenever
        storage capacity allows, particularly in collaborative projects or
        long-term archival workflows where later validation of CTF quality may
        become necessary.

        Final Perspective

        Accurate transfer of CTF metadata is an essential component of reliable
        cryo-EM data interoperability. By preserving the relationship between
        microscope imaging conditions and downstream reconstruction workflows,
        this protocol helps ensure that exported datasets remain scientifically
        consistent, reproducible, and ready for continued analysis in RELION.
    """

    _label = 'export ctf'
    _devStatus = PROD
    CTF_STAR_FILE = 'micrographs_ctf_%06d.star'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        
        form.addSection(label='Input')

        form.addParam('inputCTF', params.PointerParam,
                      pointerClass="SetOfCTF",
                      label='Input CTF',
                      help='Select set of CTF that you want to export.')

        form.addParam('micrographSource', params.EnumParam,
                      choices=['used for CTF estimation', 'other'],
                      default=0, important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Micrographs source',
                      help='By default the micrograph used to create the'
                           'exported STAR files are those used for the CTF '
                           'estimation. You can selected *other* to use a '
                           'different set of micrographs (e.g., dose weighted)')

        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      condition='micrographSource == 1',
                      important=True, label='Input micrographs',
                      help='Select the SetOfMicrographs from which to extract.')

        form.addParallelSection(threads=0, mpi=0)
            
    # -------------------------- INSERT steps functions -----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(self.writeCtfStarStep, needsGPU=False)
        
    def writeCtfStarStep(self):
        pwutils.cleanPath(self._getExportPath())
        pwutils.makePath(self._getExportPath())
        inputCTF = self.inputCTF.get()

        if self.micrographSource == 0:  # same as CTF estimation
            ctfMicSet = inputCTF.getMicrographs()
        else:
            ctfMicSet = self.inputMicrographs.get()

        micSet = SetOfMicrographs(filename=':memory:')

        psd = inputCTF.getFirstItem().getPsdFile()
        hasPsd = psd and os.path.exists(psd)

        if hasPsd:
            psdPath = self._getExportPath('PSD')
            pwutils.makePath(psdPath)
            self.info(f"Writing PSD files to {psdPath}")

        for ctf in inputCTF:
            # Get the corresponding micrograph
            mic = ctfMicSet[ctf.getObjId()]
            if mic is None:
                self.warning(f"Skipping CTF id: {ctf.getObjId()}, it is missing from input "
                             f"micrographs. ")
                continue

            micFn = mic.getFileName()
            if not os.path.exists(micFn):
                self.warning(f"Skipping micrograph {micFn}, it does not exists.")
                continue

            mic2 = mic.clone()
            mic2.setCTF(ctf)
            if hasPsd:
                psdFile = ctf.getPsdFile()
                newPsdFile = os.path.join(psdPath,
                                          '%s_psd.mrc' % pwutils.removeExt(mic.getMicName()))
                if not os.path.exists(psdFile):
                    self.warning(f"PSD file {psdFile} does not exits\n"
                                 f"Skipping micrograph {micFn}")
                    continue
                pwutils.copyFile(psdFile, newPsdFile)
                # PSD path is relative to Export dir
                newPsdFile = os.path.relpath(newPsdFile, self._getExportPath())
                ctf.setPsdFile(newPsdFile)
            else:
                # remove pointer to non-existing psd file
                ctf.setPsdFile(None)
            micSet.append(mic2)

        self.info(f"Writing set: {inputCTF} to: {self._getStarFile()}")

        micDir = self._getExportPath('Micrographs')
        pwutils.makePath(micDir)
        starWriter = convert.createWriter(rootDir=self._getExportPath(),
                                          outputDir=micDir,
                                          useBaseName=True)
        starWriter.writeSetOfMicrographs(micSet, self._getStarFile())

    # -------------------------- INFO functions -------------------------------

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

    def _getExportPath(self, *paths):
        return os.path.join(self._getPath('Export'), *paths)

    def _getStarFile(self):
        return self._getExportPath(self.CTF_STAR_FILE % self.getObjId())

    def _getPixelSize(self):
        if self.micrographSource == 0:  # same as CTF estimation
            ctfMicSet = self.inputCTF.get().getMicrographs()
        else:
            ctfMicSet = self.inputMicrographs.get()

        return ctfMicSet.getSamplingRate()
