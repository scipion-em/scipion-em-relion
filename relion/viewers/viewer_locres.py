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

from pyworkflow.viewer import ProtocolViewer
from pwem.emlib.image import ImageHandler

from .viewer_base import *
from ..protocols import ProtRelionLocalRes

binaryCondition = ('(colorMap == %d) ' % COLOR_OTHER)


class RelionLocalResViewer(ProtocolViewer):
    """ Visualization of Relion local resolution results. """

    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtRelionLocalRes]
    _label = 'viewer localres'

    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('Colored resolution')
        group.addParam('colorMap', params.EnumParam,
                       choices=list(COLOR_CHOICES.values()),
                       default=COLOR_JET,
                       label='Color map',
                       help='Select the color map to apply to the resolution '
                            'map. http://matplotlib.org/1.3.0/examples/color/'
                            'colormaps_reference.html.')
        group.addParam('otherColorMap', params.StringParam, default='jet',
                       condition=binaryCondition,
                       label='Customized Color map',
                       help='Name of a color map to apply to the resolution '
                            'map. Valid names can be found at http://'
                            'matplotlib.org/1.3.0/examples/color/'
                            'colormaps_reference.html')

        group = form.addGroup('Slices')
        group.addParam('sliceAxis', params.EnumParam, default=AX_Z,
                       choices=['x', 'y', 'z'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Slice axis')
        group.addParam('doShowVolumeSlices', params.LabelParam,
                       label="Show volume slices")
        group.addParam('doShowChimera', params.LabelParam,
                       label="Show colored map in Chimera", default=True)

    def _getVisualizeDict(self):
        self.protocol._createFilenameTemplates()
        return {
            'doShowVolumeSlices': self._showVolumeSlices,
            'doShowChimera': self._showChimera,
        }

    # =============================================================================
    # doShowVolumeSlices
    # =============================================================================
    def _showVolumeSlices(self, param=None):
        imageFile = self.protocol._getFileName('resolMap')
        imgData, minRes, maxRes = self._getImgData(imageFile)

        xplotter = RelionPlotter(x=2, y=2, mainTitle="Local Resolution Slices "
                                                     "along %s-axis."
                                                     % self._getAxis())
        for i in range(4):
            slice = self._getSlice(i + 1, imgData)
            a = xplotter.createSubPlot("Slice %s" % slice, '', '')
            matrix = self._getSliceImage(imgData, i + 1, self._getAxis())
            plot = xplotter.plotMatrix(a, matrix, minRes, maxRes,
                                       cmap=self._getColorName(),
                                       interpolation="nearest")
        xplotter.getColorBar(plot)
        return [xplotter]

    # =============================================================================
    # showChimera
    # =============================================================================
    def _showChimera(self, param=None):
        cmdFile = self.protocol._getExtraPath('chimera_local_res.cmd')
        self._createChimeraScript(cmdFile)
        view = ChimeraView(cmdFile)
        return [view]

    # =============================================================================
    # Utils Functions
    # =============================================================================
    def _getAxis(self):
        return self.getEnumText('sliceAxis')

    def _getImgData(self, imgFile):
        import numpy as np
        img = ImageHandler().read(imgFile + ":mrc")
        imgData = img.getData()

        maxRes = np.amax(imgData)
        imgData2 = np.ma.masked_where(imgData < 0.1, imgData, copy=True)
        minRes = np.amin(imgData2)

        return imgData2, minRes, maxRes

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

    def _getColorName(self):
        if self.colorMap.get() != COLOR_OTHER:
            return COLOR_CHOICES[self.colorMap.get()]
        else:
            return self.otherColorMap.get()

    def _createChimeraScript(self, scriptFile):
        import pyworkflow.gui.plotter as plotter
        fhCmd = open(scriptFile, 'w')
        imageFile = os.path.abspath(self.protocol._getFileName('resolMap'))

        _, minRes, maxRes = self._getImgData(imageFile)

        stepColors = self._getStepColors(minRes, maxRes)
        colorList = plotter.getHexColorList(stepColors, self._getColorName())

        fnVol = os.path.abspath(self.protocol._getFileName('outputVolume'))

        fhCmd.write("background solid white\n")

        fhCmd.write("open %s\n" % fnVol)
        fhCmd.write("open %s\n" % imageFile)

        sampRate = self.protocol.outputVolume.getSamplingRate()
        fhCmd.write("volume #0 voxelSize %s\n" % (str(sampRate)))
        fhCmd.write("volume #1 voxelSize %s\n" % (str(sampRate)))
        fhCmd.write("volume #1 hide\n")

        scolorStr = ''
        for step, color in zip(stepColors, colorList):
            scolorStr += '%s,%s:' % (step, color)
        scolorStr = scolorStr[:-1]
        line = ("scolor #0 volume #1 perPixel false cmap " + scolorStr + "\n")
        fhCmd.write(line)

        scolorStr2 = ''
        for step, color in zip(stepColors, colorList):
            indx = stepColors.index(step)
            if (indx % 4) != 0:
                scolorStr2 += '" " %s ' % color
            else:
                scolorStr2 += '%s %s ' % (step, color)
        line = ("colorkey 0.01,0.05 0.02,0.95 labelColor None "
                + scolorStr2 + " \n")
        fhCmd.write(line)
        fhCmd.close()

    def _getStepColors(self, minRes, maxRes, numberOfColors=13):
        inter = (maxRes - minRes) / (numberOfColors - 1)
        rangeList = []
        for step in range(0, numberOfColors):
            rangeList.append(round(minRes + step * inter, 2))
        return rangeList
