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
    """ Export a SetOfCTF to a Relion STAR file. """

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
        self._insertFunctionStep('writeCtfStarStep')
        
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
            print("Writing PSD files to %s" % psdPath)

        for ctf in inputCTF:
            # Get the corresponding micrograph
            mic = ctfMicSet[ctf.getObjId()]
            if mic is None:
                print("Skipping CTF id: %s, it is missing from input "
                      "micrographs. " % ctf.getObjId())
                continue

            micFn = mic.getFileName()
            if not os.path.exists(micFn):
                print("Skipping micrograph %s, it does not exists. " % micFn)
                continue

            mic2 = mic.clone()
            mic2.setCTF(ctf)
            if hasPsd:
                psdFile = ctf.getPsdFile()
                newPsdFile = os.path.join(psdPath,
                                          '%s_psd.mrc' % pwutils.removeExt(mic.getMicName()))
                if not os.path.exists(psdFile):
                    print("PSD file %s does not exits" % psdFile)
                    print("Skipping micrograph %s" % micFn)
                    continue
                pwutils.copyFile(psdFile, newPsdFile)
                # PSD path is relative to Export dir
                newPsdFile = os.path.relpath(newPsdFile, self._getExportPath())
                ctf.setPsdFile(newPsdFile)
            else:
                # remove pointer to non-existing psd file
                ctf.setPsdFile(None)
            micSet.append(mic2)

        print("Writing set: %s to: %s" % (inputCTF, self._getStarFile()))

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
