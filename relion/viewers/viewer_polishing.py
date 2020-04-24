# ******************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
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
from ..protocols import ProtRelionBayesianPolishing


class RelionPolishViewer(ProtocolViewer):
    """ Visualization of Relion bayesian polishing results. """
    _targets = [ProtRelionBayesianPolishing]
    _environments = [DESKTOP_TKINTER]

    _label = 'viewer polishing'

    def setProtocol(self, protocol):
        ProtocolViewer.setProtocol(self, protocol)
        self.__defineParams(self._form)
        self._createVarsFromDefinition()

    def _defineParams(self, form):
        self._env = os.environ.copy()
        self._form = form

    def __defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('3D analysis')

        group.addParam('showParticles', params.LabelParam,
                       default=True, label='Display shiny particles',
                       help='')
        group.addParam('guinierPlot', params.LabelParam,
                       default=True, label='Display polishing scale-factors')
        group.addParam('bfactorPlot', params.LabelParam,
                       default=True, label='Display polishing B-factors')

    def _getVisualizeDict(self):
        self._load()
        return {'showParticles': self._showParticles,
                'guinierPlot': lambda paramName: self._showBFactorPlot(key='guinier'),
                'bfactorPlot': lambda paramName: self._showBFactorPlot(key='bfactor'),
                }

    def _showParticles(self, paramName=None):
        views = []
        fn = self.protocol._getFileName('shiny')
        if not pwutils.exists(fn):
            raise Exception("Missing data star file '%s'")
        v = self.createScipionPartView(fn)
        views.append(v)

        return views

    def _showBFactorPlot(self, key):
        if key == 'guinier':
            label = 'rlnFittedInterceptGuinierPlot'
        else:  # bfactor
            label = 'rlnBfactorUsedForSharpening'
        gridsize = [1, 1]
        md.activateMathExtensions()

        xplotter = RelionPlotter(x=gridsize[0], y=gridsize[1],
                                 windowTitle='Polishing scale-factors')
        a = xplotter.createSubPlot("", 'Movie frame number',
                                   label,
                                   yformat=False)
        legends = []
        modelStar = self.protocol._getFileName('bfactors')
        if pwutils.exists(modelStar):
            self._plotBfactors(a, modelStar, label)
            legends.append(label)

        xplotter.showLegend(legends)
        a.grid(True)

        return [xplotter]

    def _plotBfactors(self, a, model, label):
        table = Table(fileName=model, tableName='perframe_bfactors')
        frame = table.getColumnValues('rlnMovieFrameNumber')
        bfactor = table.getColumnValues(label)
        a.plot(frame, bfactor)

    def _load(self):
        self.protocol._initialize()  # Load filename templates

    def createScipionPartView(self, filename):
        inputParticlesId = self.protocol.inputParticles.get().strId()

        labels = 'enabled id _size _filename _transform._matrix'
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels, showj.RENDER: '_filename',
                      'labels': 'id',
                      }
        return ObjectView(self._project,
                          self.protocol.strId(), filename, other=inputParticlesId,
                          env=self._env,
                          viewParams=viewParams)
