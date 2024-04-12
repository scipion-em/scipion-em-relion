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
from pyworkflow.gui.plotter import Plotter

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
                       default=True, label='Display shiny particles')
        group.addParam('guinierPlot', params.LabelParam,
                       default=True, label='Display polishing scale-factors')
        group.addParam('bfactorPlot', params.LabelParam,
                       default=True, label='Display polishing B-factors')

    def _getVisualizeDict(self):
        self._load()
        return {'showParticles': self._showParticles,
                'guinierPlot': lambda paramName: self._showBFactorPlot(key='guinier'),
                'bfactorPlot': lambda paramName: self._showBFactorPlot(key='bfactor')
                }

    def _showParticles(self, paramName=None):
        views = []
        if getattr(self.protocol, 'outputParticles', None) is not None:
            fn = self.protocol.outputParticles.getFileName()
            v = self.createScipionPartView(fn)
            views.append(v)
        return views

    def _showBFactorPlot(self, key):
        if key == 'guinier':
            label = 'rlnFittedInterceptGuinierPlot'
            title = "Polishing scale-factors"
        else:  # bfactor
            label = 'rlnBfactorUsedForSharpening'
            title = "Polishing B-factors"

        modelStar = self.protocol._getFileName('bfactors')
        if os.path.exists(modelStar):
            table = Table(fileName=modelStar, tableName='perframe_bfactors')
            frame = table.getColumnValues('rlnMovieFrameNumber')
            bfactor = map(float, table.getColumnValues(label))

            plotter = Plotter()
            figure = plotter.getFigure()
            a = figure.add_subplot(111)
            a.grid(True)
            a.set_xlabel('Movie frame number')
            a.set_ylabel(label)
            a.plot(frame, list(bfactor))
            a.set_title(title)
            plotter.tightLayout()

            return [plotter]

    def _load(self):
        self.protocol._createFilenameTemplates()  # Load filename templates

    def createScipionPartView(self, filename):
        inputParticlesId = self.protocol.inputParticles.get().strId()
        labels = 'enabled id _size _filename _transform._matrix'
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels, showj.RENDER: '_filename',
                      'labels': 'id'}
        return ObjectView(self._project,
                          self.protocol.strId(), filename,
                          other=inputParticlesId,
                          env=self._env, viewParams=viewParams)
