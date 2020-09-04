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
from pwem.viewers import LocalResolutionViewer
from pwem.viewers.viewer_chimera import mapVolsWithColorkey
from pwem.wizards import ColorScaleWizardBase
from pyworkflow.viewer import ProtocolViewer

from .viewer_base import *
from ..protocols import ProtRelionLocalRes



class RelionLocalResViewer(LocalResolutionViewer):
    """ Visualization of Relion local resolution results. """

    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtRelionLocalRes]
    _label = 'viewer localres'

    def __init__(self, **kwargs):
        LocalResolutionViewer.__init__(self, **kwargs)
        self.protocol._createFilenameTemplates()

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        groupColor = form.addGroup('Colored resolution')

        group = form.addGroup('Slices')
        group.addParam('sliceAxis', params.EnumParam, default=AX_Z,
                       choices=['x', 'y', 'z'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Slice axis')
        group.addParam('doShowVolumeSlices', params.LabelParam,
                       label="Show volume slices")
        groupColor.addParam('doShowChimera', params.LabelParam,
                       label="Show colored map in Chimera", default=True)

        _, minRes, maxRes, _ = self.getImgData(self.getResolutionVolumeFileName())
        ColorScaleWizardBase.defineColorScaleParams(groupColor, defaultHighest=maxRes, defaultLowest=minRes)

    def _getVisualizeDict(self):

        return {
            'doShowVolumeSlices': self._showVolumeSlices,
            'doShowChimera': self._showChimera,
        }

    # =============================================================================
    # doShowVolumeSlices
    # =============================================================================
    def _showVolumeSlices(self, param=None):
        imageFile = self.getResolutionVolumeFileName()
        imgData, _, _, _ = self.getImgData(imageFile)

        xplotter = RelionPlotter(x=2, y=2, mainTitle="Local Resolution Slices "
                                                     "along %s-axis."
                                                     % self._getAxis())
        for i in range(4):
            slice = self._getSlice(i + 1, imgData)
            a = xplotter.createSubPlot("Slice %s" % slice, '', '')
            matrix = self._getSliceImage(imgData, i + 1, self._getAxis())
            plot = xplotter.plotMatrix(a, matrix, self.lowest.get(), self.highest.get(),
                                       cmap=self._getColorName(),
                                       interpolation="nearest")
        xplotter.getColorBar(plot)
        return [xplotter]

    def getResolutionVolumeFileName(self):
        return self.protocol._getFileName('resolMap')

    # =============================================================================
    # showChimera
    # =============================================================================
    def _showChimera(self, param=None):

        fnResVol = self.getResolutionVolumeFileName()

        vol = self.protocol.protRefine.get().outputVolume

        fnOrigMap = vol.getFileName()
        sampRate = vol.getSamplingRate()

        cmdFile = self.protocol._getExtraPath('chimera_resolution_map.py')
        self.createChimeraScript(cmdFile, fnResVol, fnOrigMap, sampRate,
                                 numColors=self.intervals.get(),
                                 lowResLimit=self.highest.get(),
                                 highResLimit=self.lowest.get())
        view = ChimeraView(cmdFile)
        return [view]

    # =============================================================================
    # Utils Functions
    # =============================================================================
    def _getAxis(self):
        return self.getEnumText('sliceAxis')

    def _getSlice(self, index, volumeData):
        return int((index + 3) * volumeData.shape[0] / 9)

    def _getSliceImage(self, volumeData, index, dataAxis):
        slice = self._getSlice(index, volumeData)
        if dataAxis == 'y':
            imgSlice = volumeData[:, slice, :]
        elif dataAxis == 'x':
            imgSlice = volumeData[:, :, slice]
        else:
            imgSlice = volumeData[slice, :, :]
        return imgSlice

