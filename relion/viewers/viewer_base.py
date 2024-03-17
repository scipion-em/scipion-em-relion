# ******************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es) [1]
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] MRC Laboratory of Molecular Biology, MRC-LMB
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

import os
from math import radians, log
from emtable import Table
import logging
logger = logging.getLogger(__name__)

import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.viewer import (DESKTOP_TKINTER, WEB_DJANGO, Viewer)
import pyworkflow.utils as pwutils
from pwem.viewers import (EmPlotter, EmProtocolViewer, showj,
                          FscViewer, DataView, ObjectView, ChimeraView,
                          ClassesView, Classes3DView, ChimeraAngDist)
from pwem.constants import ALIGN_PROJ, NO_INDEX
from pwem.objects import FSC

from relion.convert.convert_utils import relionToLocation
from ..protocols import (ProtRelionClassify2D, ProtRelionClassify3D,
                         ProtRelionRefine3D, ProtRelionInitialModel,
                         ProtRelionSelectClasses2D)
from ..constants import *


class RelionPlotter(EmPlotter):
    """ Class to create several plots. """

    def __init__(self, x=1, y=1, mainTitle="", **kwargs):
        EmPlotter.__init__(self, x, y, mainTitle, **kwargs)

    def plotMdAngularDistribution(self, title, angularMd, tableName=None, color='blue'):
        """Create an special type of subplot, representing the angular
        distribution of weight projections. A metadata should be provided containing
        labels: RLN_ORIENT_ROT, RLN_ORIENT_TILT """

        table = Table(fileName=angularMd, tableName=tableName)
        rot = radians(table.getColumnValues('rlnAngleRot'))
        tilt = radians(table.getColumnValues('rlnAngleTilt'))
        self.plotAngularDistribution(title, rot, tilt)

    def plotMd(self, mdTable, mdLabelX, mdLabelY, color='g', **args):
        """ plot metadata columns mdLabelX and mdLabelY
            if nbins is in args then and histogram over y data is made
        """
        if mdLabelX:
            xx = []
        else:
            xx = range(1, mdTable.size() + 1)
        yy = []
        if mdLabelX:
            xx += mdTable.getColumnValues(mdLabelX)
        yy += mdTable.getColumnValues(mdLabelY)

        nbins = args.pop('nbins', None)
        if nbins is None:
            self.plotData(xx, yy, color, **args)
        else:
            self.plotHist(yy, nbins, color, **args)

    def plotMdFile(self, mdFilename, mdLabelX, mdLabelY, color='g', **args):
        """ plot metadataFile columns mdLabelX and mdLabelY
            if nbins is in args then and histogram over y data is made
        """
        table = Table(fileName=mdFilename)
        self.plotMd(table, mdLabelX, mdLabelY, **args)


def protected_show(showFunc):
    def protectedShowFunc(self, paramName=None):
        try:
            return showFunc(self, paramName=paramName)
        except Exception as e:
            self._errors = [str(e)]
            return self._showErrors()

    return protectedShowFunc


class RelionViewer(EmProtocolViewer):
    """ Visualization of Relion results. """
    _targets = [ProtRelionClassify2D, ProtRelionClassify3D,
                ProtRelionRefine3D, ProtRelionInitialModel]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    _label = 'viewer'

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Visualization')
        form.addParam('viewIter', params.EnumParam,
                      choices=['last', 'selection'], default=ITER_LAST,
                      display=params.EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize",
                      help="""
*last*: only the last iteration will be visualized.
*selection*: you may specify a range of iterations.
Examples:
"1,5-8,10" -> [1,5,6,7,8,10]
"2,6,9-11" -> [2,6,9,10,11]
"2 5, 6-8" -> [2,5,6,7,8]                      
                           """)
        form.addParam('iterSelection', params.NumericRangeParam,
                      condition='viewIter==%d' % ITER_SELECTION,
                      label="Iterations list",
                      help="Write the iteration list to visualize.")

        changesLabel = 'Changes in Offset and Angles'

        group = form.addGroup('Particles')
        if self.protocol.IS_CLASSIFY:

            group.addParam('showImagesInClasses', params.LabelParam,
                           label='Show classification in Scipion', important=True,
                           help='Display each class with the number of particles assigned. \n'
                                '*Note1*: The images of one class can be shown by \n'
                                'right-click on the class and select "Open images".\n'
                                '*Note2*: This option convert the Relion star file to\n'
                                'Scipion format and can take several minutes if the \n'
                                'number of particles is high.')
            group.addParam('showClassesOnly', params.LabelParam,
                           label='Show classes only (*_model.star)',
                           help='Display the classes directly form the *_model.star file.')
            changesLabel = 'Changes in Offset, Angles and Classes'
        else:
            group.addParam('showImagesAngularAssignment', params.LabelParam,
                           label='Particles angular assignment')
        group.addParam('showOptimiserFile', params.LabelParam,
                       label='Show *_optimiser.star file')

        if self.protocol.IS_3D:
            group = form.addGroup('Volumes')

            if self._hasClasses():
                group.addParam('showClasses3D', params.EnumParam,
                               default=CLASSES_ALL,
                               choices=['all', 'selection'],
                               display=params.EnumParam.DISPLAY_HLIST,
                               label='3D Class to visualize',
                               help='')
                group.addParam('class3DSelection', params.NumericRangeParam,
                               default='1',
                               condition='showClasses3D == %d' % CLASSES_SEL,
                               label='Classes list',
                               help='')
            else:
                if self.protocol.IS_3D_INIT:
                    group.addHidden('showHalves', params.IntParam, default=3)
                else:
                    group.addParam('showHalves', params.EnumParam, default=0,
                                   choices=['half1', 'half2', 'both', 'final'],
                                   label='Volume to visualize',
                                   help='Select which half do you want to visualize.')

            group.addParam('displayVol', params.EnumParam,
                           choices=['slices', 'chimera'],
                           default=VOLUME_SLICES,
                           display=params.EnumParam.DISPLAY_HLIST,
                           label='Display volume with',
                           help='*slices*: display volumes as 2D slices along z axis.\n'
                                '*chimera*: display volumes as surface with Chimera.')
            group.addParam('displayAngDist', params.EnumParam,
                           choices=['2D plot', 'chimera'],
                           default=ANGDIST_2DPLOT,
                           display=params.EnumParam.DISPLAY_HLIST,
                           label='Display angular distribution',
                           help='*2D plot*: display angular distribution as interactive 2D in matplotlib.\n'
                                '*chimera*: display angular distribution using Chimera with red spheres.')
            group.addParam('spheresScale', params.IntParam, default=-1,
                           expertLevel=LEVEL_ADVANCED,
                           condition='displayAngDist == %d' % ANGDIST_CHIMERA,
                           label='Spheres distance',
                           help='If the value is -1 then the distance is set '
                                'to 0.75 * xVolDim')

            group = form.addGroup('Resolution')
            group.addHidden('figure', params.EnumParam, default=0,
                            choices=['new', 'active'])
            group.addParam('resolutionPlotsSSNR', params.LabelParam,
                           default=True,
                           label='Display SSNR plots',
                           help='Display signal to noise ratio plots (SSNR)')
            if not self.protocol.IS_CLASSIFY and not self.protocol.IS_3D_INIT:
                group.addParam('resolutionPlotsFSC', params.LabelParam,
                               default=True,
                               label='Display resolution plots (FSC)',
                               help='')
                group.addParam('resolutionThresholdFSC', params.FloatParam,
                               default=0.143,
                               expertLevel=LEVEL_ADVANCED,
                               label='Threshold in resolution plots',
                               help='')

        form.addSection('Overall')
        form.addParam('showPMax', params.LabelParam, default=True,
                      label="Show average PMax",
                      help='Average (per class) of the maximum value\n '
                           'of normalized probability function')
        form.addParam('showChanges', params.LabelParam, default=True,
                      label=changesLabel,
                      help='Visualize changes in orientation, offset and\n '
                           'number images assigned to each class')
        if self.protocol.IS_CLASSIFY:
            form.addParam('plotClassDistribution', params.LabelParam,
                          default=True,
                          label='Plot class distribution over iterations',
                          help='Plot each class distribution over iterations as '
                               'bar plots.')

    def _getVisualizeDict(self):
        self._load()
        visualizeDict = {
            'showImagesInClasses': self._showImagesInClasses,
            'showClassesOnly': self._showClassesOnly,
            'showImagesAngularAssignment': self._showImagesAngularAssignment,
            'showOptimiserFile': self._showOptimiserFile,
            'showLL': self._showLL,
            'showPMax': self._showPMax,
            'showChanges': self._showChanges,
            'displayVol': self._showVolumes,
            'displayAngDist': self._showAngularDistribution,
            'resolutionPlotsSSNR': self._showSSNR,
            'resolutionPlotsFSC': self._showFSC,
            'plotClassDistribution': self._plotClassDistribution,
        }

        # If the is some error during the load, just show that instead
        # of any viewer
        if self._errors:
            for k in visualizeDict.keys():
                visualizeDict[k] = self._showErrors

        return visualizeDict

    def _showErrors(self, param=None):
        views = []
        self.errorList(self._errors, views)
        return views

    def _viewAll(self, *args):
        pass

    # =============================================================================
    # showImagesInClasses
    # =============================================================================
    def _getZoom(self):
        # Ensure that classes are shown at least at 128 px to
        # properly see the rlnClassDistribution label.
        dim = self.protocol.inputParticles.get().getDim()[0]
        if dim < 128:
            zoom = 128 * 100 // dim
        else:
            zoom = 100
        return zoom

    def _showImagesInClasses(self, paramName=None):
        """ Read Relion _data.star images file and
        generate a new metadata with the Xmipp classification standard:
        a 'classes' block and a 'class00000?_images' block per class.
        If the new metadata was already written, it is just shown.
        """
        views = []
        if (self.viewIter == ITER_LAST and
                getattr(self.protocol, 'outputClasses', None) is not None):
            fn = self.protocol.outputClasses.getFileName()
            v = self.createScipionView(fn)
            views.append(v)
        else:
            for it in self._iterations:
                fn = self.protocol._getIterClasses(it)
                v = self.createScipionView(fn)
                views.append(v)

        return views

    def _showClassesOnly(self, paramName=None):
        views = []
        viewParams = {showj.MODE: showj.MODE_GALLERY,
                      showj.RENDER: 'rlnReferenceImage',
                      showj.SORT_BY: 'rlnClassDistribution desc',
                      showj.LABELS: 'rlnClassDistribution',
                      showj.ZOOM: str(self._getZoom())
                      }

        for it in self._iterations:
            modelFile = self.protocol._getFileName('model', iter=it)
            v = self.createDataView('model_classes@' + modelFile,
                                    viewParams=viewParams)
            views.append(v)
        return views

    # =============================================================================
    # showImagesAngularAssignment
    # =============================================================================
    @protected_show
    def _showImagesAngularAssignment(self, paramName=None):
        views = []

        for it in self._iterations:
            fn = self.protocol._getIterData(it, alignType=ALIGN_PROJ)
            if not os.path.exists(fn):
                raise FileNotFoundError(f"Missing data star file '{fn}'.\n"
                                        "Please select a valid iteration.")
            v = self.createScipionPartView(fn)
            views.append(v)

        return views

    @protected_show
    def _showOptimiserFile(self, paramName=None):
        views = []

        for it in self._iterations:
            optimiserFile = self.protocol._getFileName('optimiser', iter=it)
            if not os.path.exists(optimiserFile):
                raise FileNotFoundError(f"Missing optimiser file '{optimiserFile}'.\n"
                                        "Please select a valid iteration. ")
            v = self.createDataView(optimiserFile)
            views.append(v)
        return views

    # =============================================================================
    # showLLRelion
    # =============================================================================
    def _showLL(self, paramName=None):
        views = []
        for it in self._iterations:
            fn = self.protocol._getIterData(it)
            views.append(self.createScipionView(fn))

        return views

    # =============================================================================
    # ShowPMax
    # =============================================================================
    def _showPMax(self, paramName=None):
        labels = ['rlnIterationNumber', 'rlnAveragePmax',
                  'rlnLogLikelihood']
        tablePMax = Table(columns=labels)

        for it in self._getAllIters():
            if it == 1:  # skip iter1 with Pmax=1
                continue
            # always list all iterations
            prefix = self.protocol.PREFIXES[0]
            fn = self.protocol._getFileName(prefix + 'model', iter=it)
            table = Table(fileName=fn, tableName='model_general')
            row = table[0]
            tablePMax.addRow(int(it), float(row.rlnAveragePmax),
                             float(row.rlnLogLikelihood))

        fn = self.protocol._getFileName('all_avgPmax')
        with open(fn, 'w') as f:
            tablePMax.writeStar(f)

        xplotter = RelionPlotter()
        xplotter.createSubPlot("Avg PMax per Iterations", "Iterations",
                               "Avg PMax")
        xplotter.plotMd(tablePMax, 'rlnIterationNumber',
                        'rlnAveragePmax')
        xplotter.showLegend(['rlnAveragePmax'])

        return [self.createDataView(fn), xplotter]

    # =============================================================================
    # Get classes info per iteration
    # =============================================================================
    def _plotClassDistribution(self, paramName=None):
        labels = ["rlnClassDistribution", "rlnAccuracyRotations",
                  "rlnAccuracyTranslationsAngst"]

        classInfo = {}
        iterations = self._getAllIters()

        for it in iterations:
            modelStar = self.protocol._getFileName('model', iter=it)
            table = Table(fileName=modelStar, tableName='model_classes')
            for row in table:
                i, fn = relionToLocation(row.rlnReferenceImage)
                if i == NO_INDEX:  # the case for 3D classes
                    # NOTE: Since there is not an proper ID value in
                    # the classes metadata, we are assuming that class X
                    # has a filename *_classXXX.mrc (as it is in Relion)
                    # and we take the ID from there
                    index = int(fn[-7:-4])
                else:
                    index = i

                if index not in classInfo:
                    classInfo[index] = {}
                    for l in labels:
                        classInfo[index][l] = []

                for l in labels:
                    classInfo[index][l].append(float(getattr(row, l)))

        xplotter = RelionPlotter()
        xplotter.createSubPlot("Classes distribution over iterations",
                               "Iterations", "Classes Distribution")

        # Empty list for each iteration
        iters = [[]] * len(iterations)

        l = labels[0]
        for index in sorted(classInfo.keys()):
            for it, value in enumerate(classInfo[index][l]):
                iters[it].append(value)

        ax = xplotter.getLastSubPlot()

        n = len(iterations)
        ind = range(n)
        bottomValues = [0] * n
        width = 0.45  # the width of the bars: can also be len(x) sequence

        def get_cmap(N):
            import matplotlib.cm as cmx
            import matplotlib.colors as colors
            """Returns a function that maps each index in 0, 1, ... N-1 to a distinct
            RGB color."""
            color_norm = colors.Normalize(vmin=0, vmax=N)  # -1)
            scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')

            def map_index_to_rgb_color(ind):
                return scalar_map.to_rgba(ind)

            return map_index_to_rgb_color

        cmap = get_cmap(len(classInfo))

        for classId in sorted(classInfo.keys()):
            values = classInfo[classId][l]
            ax.bar(ind, values, width, label='class %s' % classId,
                   bottom=bottomValues, color=cmap(classId))
            bottomValues = [a + b for a, b in zip(bottomValues, values)]

        ax.get_xaxis().set_ticks([i + 0.25 for i in ind])
        ax.get_xaxis().set_ticklabels([str(i) for i in ind])
        ax.legend(loc='upper left', fontsize='xx-small')

        return [xplotter]

    # =============================================================================
    # ShowChanges
    # =============================================================================
    def _showChanges(self, paramName=None):
        labels = ['rlnIterationNumber'] + self.protocol.CHANGE_LABELS
        tableChanges = Table(columns=labels)

        logger.info("Computing average changes in offset, angles, and class membership")
        for it in self._getAllIters():
            fn = self.protocol._getFileName('optimiser', iter=it)
            if not os.path.exists(fn):
                continue
            logger.info(f"Computing data for iteration {it:03d}")
            fn = self.protocol._getFileName('optimiser', iter=it)
            table = Table(fileName=fn, tableName='optimiser_general', types=LABELS_DICT)
            row = table[0]
            cols = [getattr(row, value) for value in self.protocol.CHANGE_LABELS]
            tableChanges.addRow(it, *cols)

        fn = self.protocol._getFileName('all_changes')

        with open(fn, 'w') as f:
            tableChanges.writeStar(f)

        return [self.createDataView(fn)]

    # =============================================================================
    # ShowVolumes
    # =============================================================================
    @protected_show
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        elif self.displayVol == VOLUME_SLICES:
            return self._createVolumesSqlite()

    def _createVolumesSqlite(self):
        """ Write an sqlite with all volumes selected for visualization. """
        path = self.protocol._getExtraPath('relion_viewer_volumes.sqlite')
        samplingRate = self.protocol.inputParticles.get().getSamplingRate()

        files = []
        volumes = self._getVolumeNames()
        for volFn in volumes:
            if not os.path.exists(volFn):
                raise FileNotFoundError(f"Missing volume file: {volFn}\n Please select "
                                        "a valid class or iteration number.")
            logger.debug(f"Adding vol: {volFn}")
            if not volFn.endswith(":mrc"):
                files.append(volFn + ":mrc")

        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project, self.protocol.strId(), path)]

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()
        cmdFile = self.protocol._getExtraPath('chimera_volumes.cxc')
        with open(cmdFile, 'w+') as f:
            for vol in volumes:
                # We assume that the chimera script will be generated
                # at the same folder as relion volumes
                localVol = os.path.basename(vol)
                vol = vol.replace(":mrc", "")
                if os.path.exists(vol):
                    f.write("open %s\n" % localVol)
            f.write('tile\n')
        view = ChimeraView(cmdFile)
        return [view]

    # =============================================================================
    # showAngularDistribution
    # =============================================================================
    @protected_show
    def _showAngularDistribution(self, paramName=None):
        views = []

        if self.displayAngDist == ANGDIST_CHIMERA:
            for it in self._iterations:
                views.append(self._createAngDistChimera(it))

        elif self.displayAngDist == ANGDIST_2DPLOT:
            for it in self._iterations:
                plot = self._createAngDist2D(it)
                if isinstance(plot, RelionPlotter):
                    views.append(plot)
        return views

    def _createAngDistChimera(self, it):
        # Common variables to use
        nparts = self.protocol.inputParticles.get().getSize()
        radius = self.spheresScale.get()
        prefixes = self._getPrefixes()

        if len(self._refsList) != 1:
            return self.infoMessage("Please select only one class to display "
                                    "angular distribution", "Input selection")
        # If just one reference we can show the angular distribution
        ref3d = self._refsList[0]
        volFn = self._getVolumeNames()[0]
        if not os.path.exists(volFn):
            raise FileNotFoundError("This class is empty. Please try with another class")

        for prefix in prefixes:
            sqliteFn = self.protocol._getFileName('projections',
                                                  iter=it, ref3d=ref3d,
                                                  half=prefix)
            if not os.path.exists(sqliteFn):
                mdOut = self._getMdOut(it, prefix, ref3d)
                self.createAngDistributionSqlite(
                    sqliteFn, nparts,
                    itemDataIterator=self._iterAngles(mdOut))
            if hasattr(self.protocol, 'outputVolumes'):
                vol = self.protocol.outputVolumes.getFirstItem()
            else:
                vol = self.protocol.outputVolume
            volOrigin = vol.getOrigin(force=True).getShifts()
            samplingRate = vol.getSamplingRate()
            return ChimeraAngDist(volFn, self.protocol._getPath(),
                                  voxelSize=samplingRate,
                                  volOrigin=volOrigin,
                                  angularDistFile=sqliteFn,
                                  spheresDistance=radius)

    def _createAngDist2D(self, it):
        # Common variables to use
        nparts = self.protocol.inputParticles.get().getSize()
        prefixes = self._getPrefixes()
        nrefs = len(self._refsList)
        n = nrefs * len(prefixes)
        gridsize = self._getGridSize(n)
        if prefixes[0] == "final":
            title = "Final"
        else:
            title = 'Iteration %d' % it

        plotter = RelionPlotter(x=gridsize[0], y=gridsize[1],
                                mainTitle=title, windowTitle="Angular Distribution")
        for prefix in prefixes:
            for ref3d in self._refsList:
                randomSet = self._getRandomSet(prefix)
                if randomSet > 0:
                    title = '%s class %d' % (prefix, ref3d)
                else:
                    title = 'class %d' % ref3d
                sqliteFn = self.protocol._getFileName('projections',
                                                      iter=it, ref3d=ref3d, half=prefix)
                if not os.path.exists(sqliteFn):
                    self.createAngDistributionSqlite(sqliteFn, nparts,
                                                     itemDataIterator=self._iterAngles(
                                                         self._getMdOut(it, prefix, ref3d)))
                plotter.plotAngularDistributionFromMd(sqliteFn, title)

        for prefix in prefixes:
            dataStar = self._getDataStar(prefix, it)
            if os.path.exists(dataStar):
                return plotter
            else:
                return

    # =============================================================================
    # plotSSNR
    # =============================================================================
    def _showSSNR(self, paramName=None):
        prefixes = self._getPrefixes()
        nrefs = len(self._refsList)
        n = nrefs * len(prefixes)
        gridsize = self._getGridSize(n)
        xplotter = RelionPlotter(x=gridsize[0], y=gridsize[1])

        for prefix in prefixes:
            for ref3d in self._refsList:
                plot_title = 'Resolution SSNR %s, for Class %s' % (prefix, ref3d)
                a = xplotter.createSubPlot(plot_title, 'Angstroms^-1', 'log(SSNR)')
                blockName = 'model_class_%d' % ref3d
                for it in self._iterations:
                    fn = self._getModelStar(prefix, it)
                    if os.path.exists(fn):
                        self._plotSSNR(a, fn, blockName, 'iter %d' % it)
                xplotter.legend()
                a.grid(True)

        return [xplotter]

    def _plotSSNR(self, a, fn, table, label):
        table = Table(fileName=fn, tableName=table)
        ssnr = map(float, table.getColumnValues('rlnSsnrMap'))
        resolution_inv = map(float, table.getColumnValues('rlnResolution'))
        ssnrDict = {k: v for (k, v) in zip(ssnr, resolution_inv)}
        ssnrNewDict = {}

        for ssnr in ssnrDict:
            # only cross by 1 is important
            if ssnr > 0.9:
                ssnrNewDict[log(ssnr)] = ssnrDict[ssnr]

        resolution_inv = list(ssnrNewDict.values())
        frc = list(ssnrNewDict.keys())
        a.plot(resolution_inv, frc, label=label)
        a.xaxis.set_major_formatter(self._plotFormatter)

    # =============================================================================
    # plotFSC
    # =============================================================================
    def _showFSC(self, paramName=None):
        logger.info(f"Showing FSC for iterations: {self._iterations}")
        threshold = self.resolutionThresholdFSC.get()

        fscViewer = FscViewer(project=self.protocol.getProject(),
                              threshold=threshold,
                              figure=None,
                              protocol=self.protocol,
                              addButton=True)
        fscSet = self.protocol._createSetOfFSCs()
        for it in self._iterations:
            model_star = self._getModelStar('half1_', it)
            if os.path.exists(model_star):
                fsc = self._plotFSC(None, model_star, 'iter %d' % it)
                fscSet.append(fsc)

        fscViewer.visualize(fscSet)

    def _plotFSC(self, a, model_star, label, legend=None):
        if legend is None:
            legend = label
        table = Table(fileName=model_star, tableName='model_class_1')
        resolution_inv = table.getColumnValues('rlnResolution')
        frc = table.getColumnValues('rlnGoldStandardFsc')
        fsc = FSC(objLabel=legend)
        fsc.setData(resolution_inv, frc)

        return fsc

    # =============================================================================
    # Utils Functions
    # =============================================================================
    def _validate(self):
        if self.lastIter is None:
            return ['There are no completed iterations.']

    def createDataView(self, filename, viewParams={}):
        return DataView(filename, env=self._env, viewParams=viewParams)

    def createScipionView(self, filename):
        labels = 'enabled id _size _representative._filename '
        labels += '_rlnclassDistribution _rlnAccuracyRotations '
        labels += '_rlnAccuracyTranslationsAngst '
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels,
                      showj.RENDER: '_representative._filename',
                      showj.SORT_BY: '_size desc',
                      showj.ZOOM: str(self._getZoom())
                      }
        inputParticlesId = self.protocol.inputParticles.get().strId()
        ViewClass = ClassesView if self.protocol.IS_2D else Classes3DView
        view = ViewClass(self._project,
                         self.protocol.strId(), filename, other=inputParticlesId,
                         env=self._env,
                         viewParams=viewParams)

        return view

    def createScipionPartView(self, filename):
        inputParticlesId = self.protocol._getInputParticles().strId()

        labels = 'enabled id _size _filename _transform._matrix'
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels, showj.RENDER: '_filename',
                      'labels': 'id',
                      }
        return ObjectView(self._project,
                          self.protocol.strId(), filename, other=inputParticlesId,
                          env=self._env,
                          viewParams=viewParams)

    def _getRange(self, var, label):
        """ Check if the range is not empty.
        :param var: The variable to retrieve the value
        :param label: the labe used for the message string
        :return: the list with the range of values, empty
        """
        value = var.get()
        if value is None or not value.strip():
            self._errors.append('Provide %s selection.' % label)
            result = []
        else:
            result = pwutils.getListFromRangeString(value)

        return result

    def _hasClasses(self):
        p = self.protocol
        return p.IS_CLASSIFY or p.IS_3D_INIT

    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        self._refsList = [1]
        self._errors = []

        if self.protocol.IS_3D and self._hasClasses():
            if self.showClasses3D == CLASSES_ALL:
                self._refsList = range(1, self.protocol.numberOfClasses.get() + 1)
            else:
                self._refsList = self._getRange(self.class3DSelection, 'classes 3d')

        self.protocol._initialize()  # Load filename templates
        self.firstIter = self.protocol._firstIter()
        self.lastIter = self.protocol._lastIter()

        halves = getattr(self, 'showHalves', None)
        if (self.viewIter.get() == ITER_LAST or
                (halves == 3 and not self.protocol.IS_3D_INIT)):
            self._iterations = [self.lastIter]
        else:
            self._iterations = self._getRange(self.iterSelection, 'iterations')

        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq)

    def _getAllIters(self):
        """ Return all the iterations.
        For most protocols this is just the range from min to max,
        but iteration numbers in initial volume protocol are not
        contiguous.
        """
        iterations = range(self.firstIter, self.lastIter + 1)

        return [it for it in iterations
                if os.path.exists(self.protocol._getFileName('optimiser',
                                                             iter=it))]

    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1 / value
        return "1/%0.2f" % inv

    def _getGridSize(self, n=None):
        """ Figure out the layout of the plots given the number of references. """
        if n is None:
            n = len(self._refsList)

        if n == 1:
            gridsize = [1, 1]
        elif n == 2:
            gridsize = [2, 1]
        else:
            gridsize = [(n + 1) // 2, 2]

        return gridsize

    def _getPrefixes(self):
        prefixes = self.protocol.PREFIXES
        halves = getattr(self, 'showHalves', None)
        if halves:
            if halves == 0:
                prefixes = ['half1_']
            elif halves == 1:
                prefixes = ['half2_']
            elif halves == 3:
                prefixes = ['final']
        return prefixes

    def _iterAngles(self, mdOut):
        for row in mdOut:
            rot = float(row.rlnAngleRot)
            tilt = float(row.rlnAngleTilt)
            yield rot, tilt

    def _getVolumePrefixes(self):
        prefixes = self._getPrefixes()
        if prefixes[0] == 'final':
            prefixes += ['final_half1_', 'final_half2_']
        if prefixes[0] == 'final' and self.protocol.IS_3D_INIT:
            prefixes = ['finalSGD']
        return prefixes

    def _getVolumeNames(self):
        vols = []
        prefixes = self._getVolumePrefixes()
        logger.debug(f"self._iterations: {self._iterations}")
        for it in self._iterations:
            for ref3d in self._refsList:
                for prefix in prefixes:
                    volFn = self.protocol._getFileName(prefix + 'volume',
                                                       iter=it, ref3d=ref3d)
                    if os.path.exists(volFn):
                        volFn = volFn.replace(":mrc", "")
                        vols.append(volFn)
                    else:
                        raise FileNotFoundError(f"Volume {volFn} does not exists.\n"
                                                "Please select a valid iteration number.")
        return vols

    def _getMdOut(self, it, prefix, ref3d):
        randomSet = self._getRandomSet(prefix)
        dataStar = self._getDataStar(prefix, it)
        mdOut = []

        table = Table(fileName=dataStar, tableName='particles')
        for row in table:
            if 0 < randomSet < 3:
                if int(row.rlnRandomSubset) == randomSet and int(row.rlnClassNumber) == ref3d:
                    mdOut.append(row)
            else:
                if int(row.rlnClassNumber) == ref3d:
                    mdOut.append(row)

        return mdOut

    def _getDataStar(self, prefix, it):
        randomSet = self._getRandomSet(prefix)
        if randomSet > 0 or self.protocol.IS_3D_INIT:
            return self.protocol._getFileName('data', iter=it)
        else:
            return self.protocol._getFileName('dataFinal')

    def _getModelStar(self, prefix, it):
        randomSet = self._getRandomSet(prefix)
        if self.protocol.IS_3D_INIT:
            return self.protocol._getFileName('model', iter=it)
        if randomSet > 0:
            return self.protocol._getFileName(prefix + 'model', iter=it)
        else:
            return self.protocol._getFileName('modelFinal')

    def _getRandomSet(self, prefix):
        if prefix == "final":
            return 0
        elif prefix == "half1_":
            return 1
        elif prefix == "half2_":
            return 2
        else:
            return 3


class RelionSelectClassesViewer(Viewer):
    _targets = [ProtRelionSelectClasses2D]

    def _visualize(self, obj, **kwargs):
        labels = 'enabled id _size _representative._filename '
        labels += '_rlnClassDistribution _rlnPredictedClassScore '
        labels += '_rlnEstimatedResolution '
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels,
                      showj.RENDER: '_representative._filename',
                      showj.SORT_BY: '_rlnPredictedClassScore desc',
                      showj.MODE: showj.MODE_MD
                      }
        prot = self.protocol
        output = prot.outputClasses
        view = ClassesView(self.getProject(), prot.strId(), output.getFileName(),
                           viewParams=viewParams)
        return [view]
