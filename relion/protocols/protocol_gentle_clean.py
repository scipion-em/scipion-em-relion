# **************************************************************************
# *
# * Authors: Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
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
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import re
from glob import glob

from pyworkflow.protocol import Protocol
from pyworkflow.constants import PROD
from pyworkflow.protocol.params import LabelParam


class ProtRelionCleanJobs(Protocol):
    """
    Performs a gentle cleanup of completed RELION processing jobs within a
    Scipion project. The protocol is designed to reduce storage usage while
    preserving the most relevant reconstruction results and maintaining the
    integrity of the project structure.

    AI Generated:

    Clean Project (ProtRelionCleanJobs) — User Manual
        Overview

        The Clean Project protocol provides a controlled strategy for removing
        unnecessary intermediate files generated during RELION processing.
        Cryo-EM refinement, classification, polishing, and correction workflows
        can produce very large numbers of temporary and iteration-dependent
        files, many of which are no longer needed once a protocol has completed
        successfully. This protocol helps maintain manageable project sizes
        without deleting biologically relevant outputs.

        Rather than permanently deleting information, the protocol moves
        selected files into a dedicated Trash folder inside the project. This
        approach provides a safer alternative to manual cleanup because users
        can still inspect or recover moved files if needed. The protocol is
        therefore particularly useful in long-term cryo-EM projects, shared
        computing facilities, and large institutional storage environments.

        General Cleaning Philosophy

        The protocol follows a conservative cleaning philosophy intended to
        preserve the final scientific results while reducing storage occupied
        by intermediate processing stages. In iterative RELION workflows,
        intermediate iterations often consume far more disk space than the
        final outputs. In most biological analyses, only the final iteration
        contains the reconstruction state that will be interpreted,
        visualized, refined further, or deposited.

        The cleanup procedure therefore retains the most recent iteration while
        relocating older intermediate iterations. This strategy preserves the
        ability to continue biological interpretation and downstream analysis
        while minimizing storage overhead.

        Supported Workflow Types

        The protocol is intended for RELION-based workflows that generate
        extensive intermediate files. Typical examples include 2D
        classification, 3D classification, 3D refinement, initial model
        generation, multibody refinement, particle polishing, motion
        correction, and CTF refinement.

        Different workflow types produce different categories of auxiliary
        files. Some protocols mainly generate diagnostic plots and metadata,
        while others produce large iterative reconstructions or temporary
        optimization files. The cleaning strategy adapts to these differences
        automatically so that biologically meaningful outputs remain available.

        Biological and Practical Considerations

        From a biological perspective, the protocol should only be applied
        after users are confident that the relevant processing stages have
        completed successfully and that no additional inspection of earlier
        iterations is required. Intermediate iterations may still contain
        useful information when studying convergence behavior, classification
        instability, or heterogeneous particle populations.

        For publication-quality studies or difficult datasets, users may wish
        to archive important intermediate results externally before cleaning.
        Although the protocol preserves final outputs, earlier iterations can
        sometimes provide insight into refinement dynamics or alternative
        structural states.

        The protocol is particularly valuable in facilities processing many
        datasets simultaneously, where uncontrolled accumulation of temporary
        files can rapidly exhaust storage resources. In such environments,
        routine cleanup improves project portability, backup efficiency, and
        long-term maintainability.

        Safety and Recovery

        A key feature of the protocol is that files are relocated rather than
        irreversibly deleted. This provides an additional level of operational
        safety for users who may later realize that some intermediate files
        remain important for troubleshooting or reproducibility.

        Because the protocol creates a dedicated Trash directory inside the
        project, users retain the possibility of manually restoring files if
        necessary. Nevertheless, users should still exercise caution before
        running cleanup operations on actively developing projects.

        Outputs and Project Organization

        After execution, the project directory becomes significantly more
        compact and easier to navigate. The retained outputs correspond to the
        latest reconstruction state and the principal scientific products of
        each workflow. Intermediate files and auxiliary diagnostics are grouped
        within the Trash folder, improving project readability while avoiding
        permanent data loss.

        In large cryo-EM studies involving multiple rounds of refinement,
        classification, or polishing, this cleanup strategy can substantially
        reduce storage requirements without compromising downstream biological
        interpretation.

        Final Perspective

        In modern cryo-EM workflows, data management is an essential component
        of reproducible structural biology. The Clean Project protocol provides
        a practical compromise between preserving scientific reproducibility and
        maintaining sustainable storage usage. By safely relocating unnecessary
        intermediate files while keeping the final biologically relevant
        results accessible, the protocol supports efficient long-term project
        maintenance in both individual and facility-scale environments.
    """
    _label = 'clean project'
    _devStatus = PROD

    def __init__(self, **kwargs):
        Protocol.__init__(self, **kwargs)
        self.moveDict = dict()

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('desc', LabelParam,
                      label="Gentle clean procedure will move all "
                            "intermediate files from finished Relion "
                            "protocols to Trash folder. For iteration-based "
                            "jobs, only the last iteration files are kept.")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.findRelionProtsStep, needsGPU=False)
        self._insertFunctionStep(self.runCleanStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def findRelionProtsStep(self):
        project = self.getProject()
        runs = project.getRuns()
        prjPath = os.path.abspath(project.getPath())

        if os.path.isdir(os.path.join(prjPath, "Trash")):
            raise FileExistsError("Folder %s already exists!!!" %
                                  os.path.join(prjPath, "Trash"))
        else:
            os.makedirs(os.path.join(prjPath, "Trash"))

        # make a dict with finished Relion protocols
        protDict = dict()
        for prot in runs:
            protCls = prot.getClassName()
            if prot.getStatus() == "finished" and protCls.startswith("ProtRelion"):
                extraDir = prot._getExtraPath()
                if protCls not in protDict:
                    protDict[protCls] = [extraDir]
                else:
                    protDict[protCls].append(extraDir)

        if protDict:
            self.info("Found the following folders:\n")
            for k, v in protDict.items():
                self.info("%s: %s" % (k, v))
        else:
            return

        fnsTemplate = {
            'ProtRelionExtractParticles': ['micrographs_*.star'],
            'ProtRelionPostprocess': ['*.eps'],
            'ProtRelionMotioncor': ['*corrected_micrographs.star',
                                    '*.log',
                                    '*.TXT'],
            'ProtRelionBayesianPolishing': ['*_FCC_cc.mrc',
                                            '*_FCC_w0.mrc',
                                            '*_FCC_w1.mrc',
                                            '*.eps',
                                            '*_shiny.star',
                                            '*_tracks.star'],
            'ProtRelionCtfRefinement': ['*_wAcc_optics-group*.mrc',
                                        '*_xyAcc_optics-group*.mrc',
                                        '*_aberr-Axx_optics-group_*.mrc',
                                        '*_aberr-Axy_optics-group_*.mrc',
                                        '*_aberr-Ayy_optics-group_*.mrc',
                                        '*_aberr-bx_optics-group_*.mrc',
                                        '*_aberr-by_optics-group_*.mrc',
                                        '*_mag_optics-group_*.mrc',
                                        '*_fit.star', '*_fit.eps'],
        }

        baseProts = ['ProtRelionClassify2D', 'ProtRelionClassify3D',
                     'ProtRelionRefine3D', 'ProtRelionInitialModel',
                     'ProtRelionMultiBody']

        for prot in protDict:
            if prot in fnsTemplate:
                for protDir in protDict[prot]:
                    for regex in fnsTemplate[prot]:
                        self.moveDict[os.path.join(prjPath, protDir, regex)] = os.path.join(prjPath, "Trash", protDir)
            elif prot in baseProts:
                for protDir in protDict[prot]:
                    files = sorted(glob(os.path.join(prjPath, protDir, "relion_?t???_*")))
                    # Move all files except for the last iteration
                    self.debug(files)
                    result = None
                    if files:
                        s = re.search(r"_?t(\d{3})_", files[-1])
                        if s:
                            # group 1 is 3 digits iteration number
                            result = "relion_[ic]t%03d_" % int(s.group(1))
                            self.info(f"I'll keep files: {os.path.join(protDir, result)}*")
                    if result:
                        for f in files:
                            match = re.search(result, f)
                            if not match:
                                self.moveDict[os.path.join(prjPath, protDir, f)] = os.path.join(prjPath, "Trash", protDir)

    def runCleanStep(self):
        if not self.moveDict:
            self.info("Did not find any files to remove.")
            return

        self.info("Running gentle clean for finished Relion protocols..")
        for k, v in self.moveDict.items():
            try:
                os.makedirs(v, exist_ok=True)
                os.system("mv %s %s 2> /dev/null" % (k, v))
            except:
                pass
        self.info("DONE!")

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        if self.isFinished():
            prjPath = os.path.abspath(self.getProject().getPath())
            return ["Files moved to: *%s*" %
                    os.path.join(prjPath, "Trash")]
        else:
            return []
