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
from pyworkflow.constants import BETA
from pyworkflow.protocol.params import LabelParam


class ProtRelionCleanJobs(Protocol):
    """
    Run Relion gentle clean procedure for the whole project.
    """
    _label = 'clean project'
    _devStatus = BETA

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
        self._insertFunctionStep("findRelionProtsStep")
        self._insertFunctionStep("runCleanStep")
        self._insertFunctionStep('createOutputStep')

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
            'ProtRelionExtractParticles': ['../micrographs_*.star'],
            'ProtRelionPostProcess': ['*masked.mrc'],
            'ProtRelionMotionCor': ['*corrected_micrographs.star',
                                    '*.log',
                                    '*.TXT'],
            'ProtRelionBayesianPolishing': ['*_FCC_cc.mrc',
                                            '*_FCC_w0.mrc',
                                            '*_FCC_w1.mrc'],
            'ProtRelionCtfRefinement': ['*_wAcc_optics-group*.mrc',
                                        '*_xyAcc_optics-group*.mrc',
                                        '*_aberr-Axx_optics-group_*.mrc',
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
                        s = re.search("_?t(\d{3})_", files[-1])
                        if s:
                            # group 1 is 3 digits iteration number
                            result = "relion_[ic]t%03d_" % int(s.group(1))
                            self.info("I'll keep files: %s from %s" % (result, protDir))
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
