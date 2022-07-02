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

import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.objects import Volume
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtAlignVolume


class ProtRelionSymmetrizeVolume(ProtAlignVolume):
    """
    Symmetrize a volume using Relion programs:
        *relion_align_symmetry* and *relion_image_handler*.
    """
    _label = 'symmetrize volume'
    _devStatus = PROD
    _possibleOutputs = {
        'outputVolumeAligned': Volume,
        'outputVolumeSymmetrized': Volume
    }
    
    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', params.PointerParam,
                      pointerClass='Volume',
                      label="Input volume", important=True,
                      help='Select the input volume to be symmetrized. ')

        form.addParam('symmetryGroup', params.StringParam, default='c1',
                      label="Symmetry",
                      help='Select which symmetry do you want to apply. ')
        
        form.addParallelSection(threads=0, mpi=0)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep', self.inputVolume.getObjId())

    # --------------------------- STEPS functions ------------------------------
    def createOutputStep(self, volId):
        sym = self.symmetryGroup.get()

        inFn = self._getPath('input_volume.mrc')
        alignedFn = self._getPath('volume_aligned_sym%s.mrc' % sym)
        symFn = self._getPath('volume_sym%s.mrc' % sym)
        pixSize = self.inputVolume.get().getSamplingRate()

        ImageHandler().convert(self.inputVolume.get(), inFn)

        self.runJob("relion_align_symmetry",
                    "--i %s --o %s --sym %s --angpix %0.5f" % (
                        inFn, alignedFn, sym, pixSize))

        self.runJob("relion_image_handler",
                    "--i %s --o %s --sym %s" % (alignedFn, symFn, sym))

        def _defineOutputVol(name, fn):
            vol = Volume()
            vol.copyInfo(self.inputVolume.get())
            vol.setLocation(fn)
            self._defineOutputs(**{name: vol})
            self._defineTransformRelation(self.inputVolume, vol)

        _defineOutputVol('outputVolumeAligned', alignedFn)
        _defineOutputVol('outputVolumeSymmetrized', symFn)

    # -------------------------- INFO functions -------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolumeSymmetrized'):
            summary.append("Output is not ready yet.")
        else:
            summary.append("Symmetry used: *%s*" % self.symmetryGroup.get())
        return summary
