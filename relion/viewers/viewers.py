# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from io import open
from os.path import exists
from math import radians, log

import pwem
import pwem.viewers.showj as showj
from pwem.emlib.image import ImageHandler
from pyworkflow.utils import cleanPath
from pyworkflow.protocol.constants import LEVEL_ADVANCED
import pyworkflow.protocol.params as params
from pyworkflow.viewer import (Viewer, ProtocolViewer,
                               DESKTOP_TKINTER, WEB_DJANGO)
from pwem.viewers import (EmPlotter, ObjectView, ChimeraView,
                          ChimeraClientView, ClassesView,
                          Classes3DView, FscViewer, DataView,
                          EmProtocolViewer)

import relion.convert as convert
from relion.convert.metadata import Table
from relion import Plugin
from ..protocols import (
    ProtRelionClassify2D, ProtRelionClassify3D, ProtRelionRefine3D,
    ProtRelionPostprocess, ProtRelionSortParticles,
    ProtRelionInitialModel, ProtRelionLocalRes, ProtRelionMotioncor)
from ..constants import *


class RelionPlotter(EmPlotter):
    """ Class to create several plots. """
    def __init__(self, x=1, y=1, mainTitle="", **kwargs):
        EmPlotter.__init__(self, x, y, mainTitle, **kwargs)
    
    def plotMdAngularDistribution(self, title, angularMd, table=None, color='blue'):
        """Create an special type of subplot, representing the angular
        distribution of weight projections. A metadata should be provided containing
        labels: RLN_ORIENT_ROT, RLN_ORIENT_TILT """

        table = Table(fileName=angularMd, tableName=table)
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
        self.plotMd(table, mdLabelX, mdLabelY, color='g', **args)


def protected_show(showFunc):
    def protectedShowFunc(self, paramName=None):
        try:
            return showFunc(self, paramName=paramName)
        except Exception as e:
            self._errors = [str(e)]
            return self._showErrors()

    return protectedShowFunc


class RelionViewer(EmProtocolViewer):
    """ Visualization of Relion results.

    (for protocols classify 2d/3d, 3d auto-refine and initial model)
    The visualization tools follow the recommendations of Relion 2.1 tutorial:
    http://www2.mrc-lmb.cam.ac.uk/groups/scheres/relion21_tutorial.pdf
    """
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
                           help='*2D plot*: display angular distribution as interative 2D in matplotlib.\n'
                                '*chimera*: display angular distribution using Chimera with red spheres.')
            group.addParam('spheresScale', params.IntParam, default=100, 
                           expertLevel=LEVEL_ADVANCED,
                           label='Spheres size',
                           help='')

            group = form.addGroup('Resolution')
            group.addParam('figure', params.EnumParam, default=0,
                           choices=['new', 'active'],
                           label='Figure',
                           display=params.EnumParam.DISPLAY_HLIST)
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
            fn = self.protocol._getIterData(it, alignType=pwem.constants.ALIGN_PROJ)
            if not os.path.exists(fn):
                raise Exception("Missing data star file '%s'. \n"
                                "Plese select a valid iteration. "
                                % fn)
            v = self.createScipionPartView(fn)
            views.append(v)
        
        return views

    @protected_show
    def _showOptimiserFile(self,  paramName=None):
        views = []
        
        for it in self._iterations:
            optimiserFile = self.protocol._getFileName('optimiser', iter=it)
            if not os.path.exists(optimiserFile):
                raise Exception("Missing optimiser file '%s'. \n"
                                "Plese select a valid iteration. "
                                % optimiserFile)
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
            # always list all iterations
            prefix = self.protocol.PREFIXES[0]
            fn = self.protocol._getFileName(prefix + 'model', iter=it)
            table = Table(fileName=fn, tableName='model_general')
            row = table[0]
            tablePMax.addRow(it, row.rlnAveragePmax,
                             row.rlnLogLikelihood)

        fn = self.protocol._getFileName('all_avgPmax')
        with open(fn, 'w') as f:
            tablePMax.writeStar(f)

        xplotter = RelionPlotter()
        xplotter.createSubPlot("Avg PMax per Iterations", "Iterations",
                               "Avg PMax")
        xplotter.plotMd(tablePMax, 'rlnIterationNumber',
                        'rlnAveragePmax', 'g')
        xplotter.showLegend(['rlnAveragePmax'])

        return [self.createDataView(fn), xplotter]

# =============================================================================
# Get classes info per iteration
# =============================================================================
    def _plotClassDistribution(self, paramName=None):
        labels = ["rlnClassDistribution", "rlnAccuracyRotations"]
        if Plugin.IS_GT30():
            labels.append("rlnAccuracyTranslationsAngst")
        else:
            labels.append("rlnAccuracyTranslations")

        iterations = range(self.firstIter, self.lastIter + 1)
        classInfo = {}

        for it in iterations:
            modelStar = self.protocol._getFileName('model', iter=it)
            table = Table(fileName=modelStar, tableName='model_classes')
            for row in table:
                i, fn = convert.relionToLocation(row.rlnReferenceImage)
                if i == pwem.constants.NO_INDEX:  # the case for 3D classes
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

            def map_index_to_rgb_color(index):
                return scalar_map.to_rgba(index)

            return map_index_to_rgb_color

        cmap = get_cmap(len(classInfo))

        for classId in sorted(classInfo.keys()):
            values = classInfo[classId][l]
            ax.bar(ind, values, width, label='class %s' % classId,
                   bottom=bottomValues, color=cmap(classId))
            bottomValues = [a+b for a, b in zip(bottomValues, values)]

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

        print("Computing average changes in offset, angles, and class membership")
        for it in self._getAllIters():
            fn = self.protocol._getFileName('optimiser', iter=it)
            if not os.path.exists(fn):
                continue
            print("Computing data for iteration; %03d" % it)
            fn = self.protocol._getFileName('optimiser', iter=it)
            table = Table(fileName=fn, tableName='optimiser_general')
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
            if not exists(volFn.replace(':mrc', '')):
                raise Exception("Missing volume file: %s\n Please select "
                                "a valid class or iteration number."
                                % volFn)
            print("Adding vol: %s" % volFn)
            files.append(volFn)

        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project, self.protocol.strId(), path)]
    
    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()
        
        if len(volumes) > 1:
            cmdFile = self.protocol._getExtraPath('chimera_volumes.cmd')
            f = open(cmdFile, 'w+')
            for volFn in volumes:
                # We assume that the chimera script will be generated
                # at the same folder than relion volumes
                vol = volFn.replace(':mrc', '')
                localVol = os.path.basename(vol)
                if exists(vol):
                    f.write("open %s\n" % localVol)
            f.write('tile\n')
            f.close()
            view = ChimeraView(cmdFile)
        else:
            view = ChimeraClientView(volumes[0])
            
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
        if radius < 0:
            radius = self.protocol.maskDiameterA.get()/2
            
        prefixes = self._getPrefixes()

        if len(self._refsList) != 1:
            return self.infoMessage("Please select only one class to display "
                                    "angular distribution", "Input selection")
        # If just one reference we can show the angular distribution
        ref3d = self._refsList[0]
        volFn = self._getVolumeNames()[0]
        if not exists(volFn.replace(":mrc", "")):
            raise Exception("This class is Empty. Please try with other class")

        for prefix in prefixes:
            sqliteFn = self.protocol._getFileName('projections',
                                                  iter=it, ref3d=ref3d,
                                                  half=prefix)
            if not exists(sqliteFn):
                mdOut = self._getMdOut(it, prefix, ref3d)
                self.createAngDistributionSqlite(
                    sqliteFn, nparts,
                    itemDataIterator=self._iterAngles(mdOut))

            return ChimeraClientView(volFn,
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
                if not exists(sqliteFn):
                    self.createAngDistributionSqlite(sqliteFn, nparts,
                                                     itemDataIterator=self._iterAngles(self._getMdOut(it, prefix, ref3d)))
                plotter.plotAngularDistributionFromMd(sqliteFn, title)
        
        for prefix in prefixes:
            dataStar = self._getDataStar(prefix, it)
            if exists(dataStar):
                return plotter
            else:
                return
    
# =============================================================================
# plotSSNR              
# =============================================================================

    def _getFigure(self):
        return None if self.figure == 0 else 'active'

    def _showSSNR(self, paramName=None):
        prefixes = self._getPrefixes()        
        nrefs = len(self._refsList)
        n = nrefs * len(prefixes)
        gridsize = self._getGridSize(n)
        md.activateMathExtensions()
        xplotter = RelionPlotter(x=gridsize[0], y=gridsize[1],
                                 figure=self._getFigure())
        
        for prefix in prefixes:
            for ref3d in self._refsList:
                plot_title = 'Resolution SSNR %s, for Class %s' % (prefix, ref3d)
                a = xplotter.createSubPlot(plot_title,
                                           'Angstroms^-1', 'log(SSNR)', yformat=False)
                blockName = 'model_class_%d' % ref3d
                for it in self._iterations:
                    fn = self._getModelStar(prefix, it)
                    if exists(fn):
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
        print("Showing FSC for iterations: ", self._iterations)
        threshold = self.resolutionThresholdFSC.get()
        md.activateMathExtensions()
        
        fscViewer = FscViewer(project=self.protocol.getProject(),
                              threshold=threshold,
                              protocol=self.protocol,
                              figure=self._getFigure(),
                              addButton=True)
        fscSet = self.protocol._createSetOfFSCs()
        for it in self._iterations:
            model_star = self._getModelStar('half1_', it)
            if exists(model_star):
                fsc = self._plotFSC(None, model_star, 'iter %d' % it)
                fscSet.append(fsc)
        fscViewer.visualize(fscSet)
        return [fscViewer]
    
    def _plotFSC(self, a, model_star, label, legend=None):
        if legend is None:
            legend = label
        table = Table(fileName=model_star, tableName='model_class_1')
        resolution_inv = table.getColumnValues('rlnResolution')
        frc = table.getColumnValues('rlnGoldStandardFsc')
        fsc = pwem.objects.FSC(objLabel=legend)
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
        labels += '_rlnclassDistribution _rlnAccuracyRotations _rlnAccuracyTranslations'
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
            result = self._getListFromRangeString(value)

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
                self._refsList = range(1, self.protocol.numberOfClasses.get()+1)
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
        iterations = range(self.firstIter, self.lastIter+1)

        return [it for it in iterations
                if os.path.exists(self.protocol._getFileName('optimiser',
                                                             iter=it))]

    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1/value
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
            gridsize = [(n+1)/2, 2]
            
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
        print("self._iterations: ", self._iterations)
        for it in self._iterations:
            for ref3d in self._refsList:
                for prefix in prefixes:
                    volFn = self.protocol._getFileName(prefix + 'volume',
                                                       iter=it, ref3d=ref3d)
                    if os.path.exists(volFn.replace(':mrc', '')):
                        vols.append(volFn)
                    else:
                        raise Exception("Volume %s does not exists. \n"
                                        "Please select a valid iteration "
                                        "number. " % volFn)
        return vols
    
    def _getMdOut(self, it, prefix, ref3d):
        randomSet = self._getRandomSet(prefix)
        dataStar = self._getDataStar(prefix, it)
        tableName = 'particles' if Plugin.IS_GT30() else None
        mdOut = []

        table = Table(fileName=dataStar, tableName=tableName)
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


class PostprocessViewer(ProtocolViewer):
    """ Visualization of Relion postprocess results. """
    _targets = [ProtRelionPostprocess]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
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
# =============================================================================
# ShowVolumes
# =============================================================================
        
    def _showVolumeShowj(self, volPath):        
        return [DataView(volPath)]
    
    def _showVolumesChimera(self, volPath):
        """ Create a chimera script to visualize selected volumes. """
        volPath = volPath.replace(':mrc', '')
        view = ChimeraClientView(volPath)
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
    
# =============================================================================
# plotFSC            
# =============================================================================
    def _getFigure(self):
        return None if self.figure == 0 else 'active'

    def _showFSC(self, paramName=None):
        threshold = self.resolutionThresholdFSC.get()

        md.activateMathExtensions()
        fscViewer = FscViewer(project=self.protocol.getProject(),
                              threshold=threshold,
                              protocol=self.protocol,
                              figure=self._getFigure(),
                              addButton=True)
        fscSet = self.protocol._createSetOfFSCs()

        modelStar = self.protocol._getExtraPath('postprocess.star')
        for label in self._getFSCLabels():
            if exists(modelStar):
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
        frc = table.getColumnValues('rlnGoldStandardFsc')
        fsc = pwem.objects.FSC(objLabel=legend)
        fsc.setData(resolution_inv, frc)

        return fsc

# =============================================================================
# plotGuinier
# =============================================================================
    def _showGuinier(self, paramName=None):
        gridsize = [1, 1]
        md.activateMathExtensions()
        
        xplotter = RelionPlotter(x=gridsize[0], y=gridsize[1],
                                 windowTitle='Guinier Plot')
        a = xplotter.createSubPlot("", 'Angstroms^-2', 'log(Amplitude)',
                                   yformat=False)
        legends = []
        modelStar = self.protocol._getExtraPath('postprocess.star')
        for label in self._getGuinerLabels():
            if exists(modelStar):
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
    
# =============================================================================
# Utils Functions
# =============================================================================
    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq) 
        
    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1/value
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


# TODO: deprecate this class
class RelionSortViewer(Viewer):
    """ Visualization of Relion sorting results."""

    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtRelionSortParticles]

    def _visualize(self, obj, **kwargs):
        views = []

        if obj.hasAttribute('outputParticles'):  # Protocol finished
            particles = obj.outputParticles
            labels = ('id enabled _index _filename _rlnSelectParticlesZscore '
                      '_coordinate._rlnAutopickFigureOfMerit _sampling '
                      '_ctfModel._defocusU _ctfModel._defocusV '
                      '_ctfModel._defocusAngle _transform._matrix')
            sortBy = '_rlnSelectParticlesZscore asc'
            strId = particles.strId()
            fn = particles.getFileName()
            views.append(ObjectView(self._project, strId, fn,
                                    viewParams={showj.ORDER: labels,
                                                showj.VISIBLE: labels,
                                                showj.SORT_BY: sortBy,
                                                showj.RENDER: '_filename'}))

            fn = obj._getExtraPath('input_particles_sorted.star')
            mdFn = md.MetaData(fn)
            # If Zscore in output images plot Zscore particle sorting
            if mdFn.containsLabel(md.RLN_SELECT_PARTICLES_ZSCORE):
                # sort output by Z-score
                mdFn.sort(md.RLN_SELECT_PARTICLES_ZSCORE)
                xplotter = RelionPlotter(windowTitle="Zscore particles sorting")
                xplotter.createSubPlot("Particle sorting", "Particle number",
                                       "Zscore")
                xplotter.plotMd(mdFn, False,
                                mdLabelY=md.RLN_SELECT_PARTICLES_ZSCORE)
                views.append(xplotter)

        return views


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
            slice = self._getSlice(i+1, imgData)
            a = xplotter.createSubPlot("Slice %s" % slice, '', '')
            matrix = self._getSliceImage(imgData, i+1, self._getAxis())
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
        cleanPath(path)
        movieSet = pwem.objects.SetOfMovies(filename=path)
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
