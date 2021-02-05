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

import pyworkflow.utils as pwutils
from pwem.objects import SetOfMovies

from .viewer_base import *
from ..protocols import ProtRelionMotioncor


class RelionMotioncorrViewer(EmProtocolViewer):
    """ Visualization of relion motioncor results. """

    _targets = [ProtRelionMotioncor]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer motioncor'

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('Micrographs')
        group.addParam('doShowMics', params.LabelParam,
                       label="Show aligned micrographs?", default=True,
                       help="Show the output aligned micrographs.")
        group.addParam('doShowMicsDW', params.LabelParam,
                       label="Show aligned DOSE-WEIGHTED micrographs?",
                       default=True,
                       help="Show the output aligned dose-weighted "
                            "micrographs.")
        group.addParam('doShowMicsPS', params.LabelParam,
                       label="Show sum of power spectra?",
                       default=True,
                       help="Show the output sum of non-dose weighted "
                            "power spectra.")

        group = form.addGroup('Movies')
        group.addParam('doShowMovies', params.LabelParam,
                       label="Show output movies?", default=True,
                       help="Show the output movies with alignment "
                            "information.")
        group.addParam('doShowFailedMovies', params.LabelParam,
                       label="Show FAILED movies?", default=True,
                       help="Create a set of failed movies "
                            "and display it.")

        group = form.addGroup('Statistics')
        group.addParam('doShowMotion', params.LabelParam,
                       label="Plot motion per frame", default=True,
                       help="Show accumulated motion for all micrographs. "
                            "Early motion default cut-off is 4 e/A2.")

    def _getVisualizeDict(self):
        self._errors = []
        visualizeDict = {'doShowMics': self._viewMics,
                         'doShowMicsDW': self._viewMicsDW,
                         'doShowMicsPS': self._viewMicsPS,
                         'doShowMovies': self._viewMovies,
                         'doShowFailedMovies': self._viewFailed,
                         'doShowMotion': self._plotMotion,
                         }
        return visualizeDict

    def _viewMics(self, param=None):
        if getattr(self.protocol, 'outputMicrographs', None) is not None:
            return [self.micrographsView(self.getProject(),
                                         self.protocol.outputMicrographs)]
        else:
            return [self.errorMessage('No output micrographs found!',
                                      title="Visualization error")]

    def _viewMicsDW(self, param=None):
        if getattr(self.protocol, 'outputMicrographsDoseWeighted', None) is not None:
            return [self.micrographsView(self.getProject(),
                                         self.protocol.outputMicrographsDoseWeighted)]
        else:
            return [self.errorMessage('No output dose-weighted micrographs found!',
                                      title="Visualization error")]

    def _viewMicsPS(self, param=None):
        labelsPs = 'enabled id _powerSpectra._filename '
        labelsPs += '_powerSpectra._samplingRate _filename'
        viewParamsPs = {showj.MODE: showj.MODE_MD,
                        showj.ORDER: labelsPs,
                        showj.VISIBLE: labelsPs,
                        showj.RENDER: '_powerSpectra._filename'
                        }

        if getattr(self.protocol, 'savePSsum', False):
            if getattr(self.protocol, 'outputMicrographs', None) is not None:
                output = self.protocol.outputMicrographs
            else:
                output = self.protocol.outputMicrographsDoseWeighted
            return [self.objectView(output, viewParams=viewParamsPs)]
        else:
            return [self.errorMessage('No output power spectra found!',
                                      title="Visualization error")]

    def _viewMovies(self, param=None):
        labelsMovie = 'enabled id _filename _samplingRate '
        labelsMovie += '_acquisition._dosePerFrame _acquisition._doseInitial '
        self.viewParamsMovie = {showj.MODE: showj.MODE_MD,
                                showj.ORDER: labelsMovie,
                                showj.VISIBLE: labelsMovie,
                                showj.RENDER: None
                                }

        if getattr(self.protocol, 'outputMovies', None) is not None:
            output = self.protocol.outputMovies
            return [self.objectView(output, viewParams=self.viewParamsMovie)]
        else:
            return [self.errorMessage('No output movies found!',
                                      title="Visualization error")]

    def _viewFailed(self, param=None):
        self.failedList = self.protocol._readFailedList()
        if not self.failedList:
            return [self.errorMessage('No failed movies found!',
                                      title="Visualization error")]
        else:
            sqliteFn = self.protocol._getPath('movies_failed.sqlite')
            self.createFailedMoviesSqlite(sqliteFn)
            return [self.objectView(sqliteFn, viewParams=self.viewParamsMovie)]

    def createFailedMoviesSqlite(self, path):
        inputMovies = self.protocol.inputMovies.get()
        pwutils.cleanPath(path)
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

    def micrographsView(self, project, micSet, other='', **kwargs):
        """ Reimplemented from base class to add extra labels. """

        RENDER_LABELS = ['thumbnail._filename', 'psdCorr._filename',
                         'plotGlobal._filename']
        EXTRA_LABELS = ['_filename', '_rlnAccumMotionTotal',
                        '_rlnAccumMotionEarly', '_rlnAccumMotionLate']

        first = micSet.getFirstItem()

        def existingLabels(labelList):
            return ' '.join([l for l in labelList if first.hasAttributeExt(l)])

        renderLabels = existingLabels(RENDER_LABELS)
        extraLabels = existingLabels(EXTRA_LABELS)
        labels = 'id enabled %s %s' % (renderLabels, extraLabels)

        viewParams = {showj.MODE: showj.MODE_MD,
                      showj.ORDER: labels,
                      showj.VISIBLE: labels,
                      showj.ZOOM: 50
                      }

        if renderLabels:
            viewParams[showj.RENDER] = renderLabels

        inputId = micSet.getObjId() or micSet.getFileName()
        return ObjectView(project,
                          inputId, micSet.getFileName(), other,
                          viewParams, **kwargs)

    def _plotMotion(self, param=None):
        if getattr(self.protocol, 'outputMicrographs', None) is not None:
            output = self.protocol.outputMicrographs
        else:
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
                                            xlabel='Micrograph',
                                            title='Accumulated motion per frame',
                                            bins=False,
                                            orderColumn='id',
                                            orderDirection='ASC')

        return [xplotter]
