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

from pyworkflow.utils import cleanPath
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.protocol.params import LabelParam
from pwem.viewers import MicrographsView, EmProtocolViewer, EmPlotter
import pwem.viewers.showj as showj
from pwem.objects import SetOfMovies

from ..protocols import ProtRelionMotioncor


class RelionMotioncorrViewer(EmProtocolViewer):
    """ Visualization of relion motioncor results. """

    _targets = [ProtRelionMotioncor]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer motioncor'

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        if self.hasMics():
            form.addParam('doShowMics', LabelParam,
                          label="Show aligned micrographs?", default=True,
                          help="Show the output aligned micrographs.")
        if self.hasDWMics():
            form.addParam('doShowMicsDW', LabelParam,
                          label="Show aligned DOSE-WEIGHTED micrographs?",
                          default=True,
                          help="Show the output aligned dose-weighted "
                               "micrographs.")
        form.addParam('doShowMovies', LabelParam,
                      label="Show output movies?", default=True,
                      help="Show the output movies with alignment "
                           "information.")
        form.addParam('doShowFailedMovies', LabelParam,
                      label="Show FAILED movies?", default=True,
                      help="Create a set of failed movies "
                           "and display it.")
        form.addParam('doShowMotion', LabelParam,
                      label="Plot motion per frame", default=True,
                      help="Show accumulated motion for all micrographs. "
                           "Early motion default cut-off is 4 e/A2.")

    def _getVisualizeDict(self):
        self._errors = []
        visualizeDict = {'doShowMovies': self._viewParam,
                         'doShowFailedMovies': self._viewParam,
                         'doShowMotion': self._plotMotion,
                         }
        if self.hasMics():
            visualizeDict.update({'doShowMics': self._viewParam})

        if self.hasDWMics():
            visualizeDict.update({'doShowMicsDW': self._viewParam})

        return visualizeDict

    def hasMics(self):
        return hasattr(self.protocol, 'outputMicrographs')

    def hasDWMics(self):
        return hasattr(self.protocol, 'outputMicrographsDoseWeighted')

    def _viewParam(self, param=None):
        labelsDef = 'enabled id _filename _samplingRate '
        labelsDef += '_acquisition._dosePerFrame _acquisition._doseInitial '
        viewParamsDef = {showj.MODE: showj.MODE_MD,
                         showj.ORDER: labelsDef,
                         showj.VISIBLE: labelsDef,
                         showj.RENDER: None
                         }
        if param == 'doShowMics':
            return [MicrographsView(self.getProject(),
                                    self.protocol.outputMicrographs)]

        elif param == 'doShowMicsDW':
            return [MicrographsView(self.getProject(),
                                    self.protocol.outputMicrographsDoseWeighted)]

        elif param == 'doShowMovies':
            if getattr(self.protocol, 'outputMovies', None) is not None:
                output = self.protocol.outputMovies
                return [self.objectView(output, viewParams=viewParamsDef)]
            else:
                return [self.errorMessage('No output movies found!',
                                          title="Visualization error")]

        elif param == 'doShowFailedMovies':
            self.failedList = self.protocol._readFailedList()
            if not self.failedList:
                return [self.errorMessage('No failed movies found!',
                                          title="Visualization error")]
            else:
                sqliteFn = self.protocol._getPath('movies_failed.sqlite')
                self.createFailedMoviesSqlite(sqliteFn)
                return [self.objectView(sqliteFn, viewParams=viewParamsDef)]

    def createFailedMoviesSqlite(self, path):
        inputMovies = self.protocol.inputMovies.get()
        cleanPath(path)
        movieSet = SetOfMovies(filename=path)
        movieSet.copyInfo(inputMovies)
        movieSet.copyItems(inputMovies,
                           updateItemCallback=self._findFailedMovies)

        movieSet.write()
        movieSet.close()

        return movieSet

    def _findFailedMovies(self, item, row):
        if item.getObjId() not in self.failedList:
            setattr(item, "_appendItem", False)

    def _plotMotion(self, param=None):
        if self.hasDWMics():
            output = self.protocol.outputMicrographsDoseWeighted
            columns = '_rlnAccumMotionTotal _rlnAccumMotionEarly _rlnAccumMotionLate'
            xplotter = EmPlotter.createFromFile(output.getFileName(), '',
                                                plotType='Plot',
                                                columnsStr=columns,
                                                colorsStr='r g b',
                                                linesStr='- - -',
                                                markersStr='. . .',
                                                xcolumn='id',
                                                ylabel='Motion per frame (A)',
                                                xlabel='Micrograph id',
                                                title='Accumulated motion per frame',
                                                bins=False,
                                                orderColumn='id',
                                                orderDirection='ASC')
            return [xplotter]
        else:
            return [self.errorMessage('Plot is available only when dose weighting is ON',
                                      title="Visualization error")]
