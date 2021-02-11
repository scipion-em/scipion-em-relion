# ******************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es) [1]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from .viewer_base import *
from ..protocols import ProtRelionPostprocess


class PostprocessViewer(ProtocolViewer):
    """ Visualization of Relion postprocess results. """
    _targets = [ProtRelionPostprocess]
    _environments = [DESKTOP_TKINTER]

    _label = 'viewer postprocess'

    def setProtocol(self, protocol):
        ProtocolViewer.setProtocol(self, protocol)
        self.__defineParams(self._form)
        self._createVarsFromDefinition()

    def _defineParams(self, form):
        self._form = form

    def __defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('3D analysis')

        group.addParam('displayVol', params.EnumParam,
                       choices=['slices', 'chimera'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       default=VOLUME_SLICES,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')
        group.addParam('displayMaskedVol', params.EnumParam,
                       choices=['slices', 'chimera'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       default=VOLUME_SLICES,
                       label='Display masked volume with',
                       help='*slices*: display masked volume as 2D slices along z axis.\n'
                            '*chimera*: display masked volume as surface with Chimera.')
        group.addParam('figure', params.EnumParam, default=0,
                       choices=['new', 'active'],
                       label='Figure',
                       display=params.EnumParam.DISPLAY_HLIST,
                       help="Plot in a new window vs the last opened one")
        group.addParam('resolutionPlotsFSC', params.EnumParam,
                       choices=['Corrected', 'Unmasked Maps', 'Masked Maps',
                                'Phase Randomized Masked Maps', 'all'],
                       default=FSC_CORRECTED,
                       display=params.EnumParam.DISPLAY_COMBO,
                       label='Display resolution plots (FSC)')
        group.addParam('resolutionThresholdFSC',
                       params.FloatParam, default=0.143,
                       expertLevel=LEVEL_ADVANCED,
                       label='Threshold in resolution plots')
        group.addParam('guinierPlots', params.LabelParam,
                       default=True, label='Display Guinier plots',
                       help='Guinier plots are the logarithm of the '
                            'amplitudes of the individual-frame '
                            'reconstructions divided by the amplitudes of '
                            'the average reconstruction from all frames '
                            'versus the square of the resolution. Linear '
                            'fits through these Guinier plots (which may '
                            'often be performed for resolution higher '
                            'than 20 Angstroms) then yield a slope (the '
                            'B-factor) and an intercept (a '
                            'resolution-independent scale-factor) which '
                            'are used to device the resolution-dependent '
                            'weghting scheme.')

    def _getVisualizeDict(self):
        self._load()
        return {'displayVol': self._showVolume,
                'displayMaskedVol': self._showMaskedVolume,
                'guinierPlots': self._showGuinier,
                'resolutionPlotsFSC': self._showFSC
                }

    # =========================================================================
    # ShowVolumes
    # =========================================================================

    def _showVolumeShowj(self, volPath):
        return [DataView(volPath)]

    def _showVolumesChimera(self, volPath):
        """ Create a chimera script to visualize selected volumes. """
        view = ChimeraView(volPath)
        return [view]

    def _showVolume(self, paramName=None):
        volPath = self.protocol._getExtraPath('postprocess.mrc:mrc')

        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera(volPath)

        elif self.displayVol == VOLUME_SLICES:
            return self._showVolumeShowj(volPath)

    def _showMaskedVolume(self, paramName=None):
        volPath = self.protocol._getExtraPath('postprocess_masked.mrc:mrc')

        if self.displayMaskedVol == VOLUME_CHIMERA:
            return self._showVolumesChimera(volPath)

        elif self.displayMaskedVol == VOLUME_SLICES:
            return self._showVolumeShowj(volPath)

    # =========================================================================
    # plotFSC
    # =========================================================================
    def _getFigure(self):
        return None if self.figure == 0 else 'active'

    def _showFSC(self, paramName=None):
        threshold = self.resolutionThresholdFSC.get()

        fscViewer = FscViewer(project=self.protocol.getProject(),
                              threshold=threshold,
                              protocol=self.protocol,
                              figure=self._getFigure(),
                              addButton=True)
        fscSet = self.protocol._createSetOfFSCs()

        modelStar = self.protocol._getExtraPath('postprocess.star')
        for label in self._getFSCLabels():
            if os.path.exists(modelStar):
                legend = self._getLegend(label)
                fsc = self._plotFSC(None, modelStar, label, legend)
                fscSet.append(fsc)
        fscViewer.visualize(fscSet)
        return [fscViewer]

    # ROB this function is duplicated
    def _plotFSC(self, a, model_star, label, legend=None):
        if legend is None:
            legend = label
        table = Table(fileName=model_star, tableName='fsc')
        resolution_inv = table.getColumnValues('rlnResolution')
        frc = table.getColumnValues(label)
        fsc = FSC(objLabel=legend)
        fsc.setData(resolution_inv, frc)

        return fsc

    # =========================================================================
    # plotGuinier
    # =========================================================================
    def _showGuinier(self, paramName=None):
        xplotter = RelionPlotter(windowTitle='Guinier Plot')
        a = xplotter.createSubPlot("", 'Angstroms^-2', 'log(Amplitude)')
        legends = []
        modelStar = self.protocol._getExtraPath('postprocess.star')
        for label in self._getGuinerLabels():
            if os.path.exists(modelStar):
                self._plotGuinier(a, modelStar, label)
                legends.append(self._getGuinerLegend(label))

        xplotter.showLegend(legends)
        a.grid(True)

        return [xplotter]

    def _plotGuinier(self, a, model, label):
        table = Table(fileName=model, tableName='guinier')
        resolSqInv = table.getColumnValues('rlnResolutionSquared')
        logAmp = table.getColumnValues(label)

        self.maxfsc = max(logAmp)
        self.minInv = min(resolSqInv)
        self.maxInv = max(resolSqInv)

        a.plot(resolSqInv, logAmp)
        a.xaxis.set_major_formatter(self._plotFormatter)

    # =========================================================================
    # Utils Functions
    # =========================================================================
    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq)

    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1 / value
        return "1/%0.2f" % inv

    def _getFSCLabels(self):
        if self.resolutionPlotsFSC.get() == 0:
            return ['rlnFourierShellCorrelationCorrected']
        elif self.resolutionPlotsFSC.get() == 1:
            return ['rlnFourierShellCorrelationUnmaskedMaps']
        elif self.resolutionPlotsFSC.get() == 2:
            return ['rlnFourierShellCorrelationMaskedMaps']
        elif self.resolutionPlotsFSC.get() == 3:
            return ['rlnCorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps']
        else:
            return ['rlnFourierShellCorrelationCorrected',
                    'rlnFourierShellCorrelationUnmaskedMaps',
                    'rlnFourierShellCorrelationMaskedMaps',
                    'rlnCorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps']

    def _getLegend(self, label):
        if label == 'rlnFourierShellCorrelationCorrected':
            return 'Corrected'
        elif label == 'rlnFourierShellCorrelationUnmaskedMaps':
            return 'Unmasked Maps'
        elif label == 'rlnFourierShellCorrelationMaskedMaps':
            return 'Masked Maps'
        else:
            return 'Phase Randomized Masked Maps'

    def _getGuinerLabels(self):
        return ['rlnLogAmplitudesOriginal',
                'rlnLogAmplitudesWeighted',
                'rlnLogAmplitudesSharpened',
                'rlnLogAmplitudesIntercept']

    def _getGuinerLegend(self, label):
        if label == 'rlnLogAmplitudesOriginal':
            return 'log(Amplitudes) Original'
        elif label == 'rlnLogAmplitudesWeighted':
            return 'log(Amplitudes) Weighted'
        elif label == 'rlnLogAmplitudesSharpened':
            return 'log(Amplitudes) Sharpened'
        else:
            return 'log(Amplitudes) Intercept'
