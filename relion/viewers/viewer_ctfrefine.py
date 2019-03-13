# ******************************************************************************
# *
# * Authors:    Roberto Marabini       (roberto@cnb.csic.es) [1]
# *             J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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

import sys
import os
import numpy as np
import matplotlib as mpl
from itertools import izip

import pyworkflow as pw
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import LabelParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pyworkflow.em.viewers.plotter import plt, EmPlotter
from pyworkflow.em.viewers import ObjectView, DataView
import pyworkflow.em.viewers.showj as showj

from relion.protocols import ProtRelionCtfRefinement


class ProtCtfREfineViewer(ProtocolViewer):
    """ viewer of Relion cTF refine"""
    _label = 'ctf refine viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtRelionCtfRefinement]

    def __init__(self,  **kwargs):
        ProtocolViewer.__init__(self,  **kwargs)
        self._micInfoList = None
        self._currentMicIndex = 0
        self._loadAnalyzeInfo()

    def _defineParams(self, form):
        self._env = os.environ.copy()
        showBeamTilt = self.protocol.doBeamtiltEstimation.get()
        form.addSection(label="Defocus")
        form.addParam('displayDefocus', LabelParam,
                      label="Show Defocus",
                      help="Display the defocus estimation.\n "
                           "Plot defocus difference vs position in micrograph\n"
                           "You may move between micrographs by typing:\:"
                           "left/right keys (move +1/-1 micrograph)\n"
                           "up/down keys (move +10/-10 micrographs)\n"
                           "page_up/page_down keys (move +100/-100 "
                           "micrographs)\n"
                           "home/end keys (move +1000/-1000 micrographs)\n")
        form.addParam('displayPlotDEfocusStdev', LabelParam,
                      label="Show defocus variance per micrograph",
                      help="Stdev of defocus per micrographs. Micrograph"
                           "with less than 7 particles are skept.")
        form.addParam('displayBeamTilt', LabelParam,
                      label="Show BeanTilt Images",
                      condition="{}".format(showBeamTilt),
                      help="Display two images: (1) phase differences from "
                           "which this estimate was derived and\n"
                           "(2) the model fitted through it c")
        form.addParam('displayParticles', LabelParam,
                      label="display Particles",
                      help="See the particles with the new CTF and Beamtilt")

    def _getVisualizeDict(self):
        return{
            'displayDefocus': self._visualizeDefocus,
            'displayBeamTilt': self._displayBeamTilt,
            'displayParticles': self._displayParticles,
            'displayPlotDEfocusStdev': self._displayPlotDefocusStdev
        }

    def onClick(self, event):
        # If click outside plot print will fail
        # because xdata is Nonetype
        try:
            ix, iy = event.xdata, event.ydata
            print 'x = %d, y = %d' % (ix, iy)
        except:
            pass

    def _displayPlotDefocusStdev(self, e=None):
        x = [mi.micId.get() for mi in self._micInfoList]
        y = [mi.stdev.get() for mi in self._micInfoList]
        self.plotter = EmPlotter(windowTitle="Defocus STdev per Micrograph")
        self.fig = self.plotter.getFigure()
        self.ax = self.plotter.createSubPlot("Defocus STdev per Micrograph",
                                             "# Micrograph", "stdev")
        im = self.ax.scatter(x, y, s=50, marker='o')
        self.fig.canvas.mpl_connect('button_press_event', self.onClick)

        self.plotter.show()

    def _visualizeDefocus(self, e=None):
        """Show matplotlib with defocus values."""
        micInfo = self._micInfoList[self._currentMicIndex]
        # disable default binding for arrows
        # because I want to use them
        # to navigate between micrographs
        # wrap in try, except because matplotlib will raise
        # an exception if the value is not in the list
        try:
            mpl.rcParams['keymap.back'].remove('left')
            mpl.rcParams['keymap.forward'].remove('right')
        except:
            pass

        self.plotter = EmPlotter(windowTitle="CTF Refinement")

        self.fig = self.plotter.getFigure()
        self.ax = self.plotter.createSubPlot(
            self._getTitle(micInfo), "Mic-Xdim", "Mic-Ydim")#, projection='3d')

        # call self.press after pressing any key
        self.fig.canvas.mpl_connect('key_press_event', self.press)

        # scatter plots loss colormap after rotation, so
        # if the user rotates the canvas repaint.
        # I guess this is a bug in matplotlib
        self.fig.canvas.mpl_connect('draw_event', self.show)

        self.show()

    def _getTitle(self, micInfo):
        return ("use arrows or page up/down or home/end to navigate.\n"
                "Displaying Mic = %s (%d)" %
                (micInfo.micName, micInfo.micId))

    def press(self, event):
        """ Change the currently shown micrograph
        when a key is pressed (increment/decrement)
        """
        sys.stdout.flush()

        if event.key == 'q':
            quit()  # Quit the whole system?

        shiftDict = {
            'left': -1, 'right': 1,
            'up': 10, 'down': -10,
            'pageup': -100, 'pagedown': 100,
            'home': -1000, 'end': 1000
        }
        shift = shiftDict[event.key]
        # Check the new micrograph index is between first and last
        newIndex = self._currentMicIndex + shift
        self._currentMicIndex = min(max(0, newIndex), len(self._micInfoList) - 1)
        self.show()

    def show(self, event=None):
        """ Draw plot """
        micInfo = self._micInfoList[self._currentMicIndex]

        if event is None:
            # I need to clear the plot otherwise
            # ols points are not removed
            self.ax.clear()
            self.ax.margins(0.05)
            self.ax.set_title(self._getTitle(micInfo))
            newFontSize = self.plotter.plot_axis_fontsize + 2
            self.ax.set_xlabel("Mic Xdim (px)", fontsize=newFontSize)
            self.ax.set_ylabel("Mic Ydim (px)", fontsize=newFontSize)
            #self.ax.set_zlabel("Defocus Difference (A)", fontsize=newFontSize)
            self.ax.grid(True)
            # FIXME: Use mic size
            self.ax.set_xlim(0, np.max(micInfo.x))
            self.ax.set_ylim(0, np.max(micInfo.y))
            # if I do not use subplots_adjust the window shrinks
            #  after redraw
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        im = self.ax.scatter(micInfo.x, micInfo.y, #zs=micInfo.defocus,
                             c=micInfo.defocusDiff, s=100, marker='o')

        self.plotter.getColorBar(im)
        self.plotter.show()

    def _displayBeamTilt(self, paramName=None):
        phaseDifferenceFn = self.protocol.fileWithPhaseDifferenceName()
        modelFitFn = self.protocol.fileWithModelFitterName()

        # TODO: how can I change the window title?
        return [DataView(phaseDifferenceFn), DataView(modelFitFn)]

    def createScipionPartView(self, filename, viewParams={}):
        inputParticlesId = self.protocol.inputParticles.get().strId()
        labels = 'enabled id _size _filename '
        labels += ' _ctfModel._defocusU _ctfModel._defocusV '

        if self.protocol.doBeamtiltEstimation:
            labels += ' _rlnBeamTiltX _rlnBeamTiltY'

        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels,
                      showj.MODE: showj.MODE_MD,
                      showj.RENDER: '_filename',
                      'labels': 'id',
                      }

        return ObjectView(self._project,
                          self.protocol.strId(),
                          filename,
                          other=inputParticlesId,
                          env=self._env,
                          viewParams=viewParams)

    def _displayParticles(self, paramName=None):
        views = []
        fn = self.protocol.outputParticles.getFileName()
        v = self.createScipionPartView(fn)
        views.append(v)
        return views

    # ------------------- UTILS functions -------------------------
    def _loadAnalyzeInfo(self):
        # Only load once
        if self._micInfoList is None:
            ctfInfoFn = self.protocol._getExtraPath('ctf_analyze.sqlite')
            if not os.path.exists(ctfInfoFn):
                ctfInfo = self._createCtfInfo(ctfInfoFn)
            else:
                ctfInfo = AnalizeCtfInfo(ctfInfoFn)
            self._micInfoList = [mi.clone() for mi in ctfInfo]
            ctfInfo.close()

    def _createCtfInfo(self, ctfInfoFn):
        ctfInfo = AnalizeCtfInfo(filename=ctfInfoFn)
        print("Generating CTF statistics...")
        ctfInfo.loadFromParticles(self.protocol.inputParticles.get(),
                                  self.protocol.outputParticles)
        return ctfInfo


class AnalizeCtfInfo:
    """ Simple class to store visualization information related
    to Micrographs and Particles CTF information after
    Relion - ctf refinement protocol.
    """
    def __init__(self, filename):
        # The classes dict needs to be updated to register local objects
        classesDict = dict(pwobj.__dict__)
        classesDict['MicInfo'] = MicInfo
        self._infoSet = pwobj.Set(filename, classesDict=classesDict)

    def addMicInfo(self, micId, x, y, defocus, defocusDiff):
        pass

    def loadFromParticles(self, inputParts, outputParts):
        # compute difference in defocus and save in database
        micInfo = None
        lastMicId = None
        infoList = []

        def _avgDefocus(p):
            dU, dV, _ = p.getCTF().getDefocus()
            return (dU + dV) / 2.0

        for p1, p2 in izip(inputParts.iterItems(orderBy=['_micId', 'id']),
                           outputParts.iterItems(orderBy=['_micId', 'id'])):
            coord = p1.getCoordinate()
            micId = coord.getMicId()

            if micId != lastMicId:
                micInfo = MicInfo(micId=micId, micName=coord.getMicName())
                infoList.append(micInfo)
                lastMicId = micId

            p1D = _avgDefocus(p1)
            p2D = _avgDefocus(p2)
            x, y = coord.getPosition()
            micInfo.addEntry(x, y, p2D, p2D - p1D)

        for micInfo in infoList:
            micInfo.computeStats()
            self._infoSet.append(micInfo)

        self._infoSet.write()

    def __iter__(self):
        for micInfo in self._infoSet:
            yield micInfo

    def close(self):
        self._infoSet.close()


class MicInfo(pwobj.OrderedObject):
    def __init__(self, **kwargs):
        pwobj.OrderedObject.__init__(self, **kwargs)

        self.micId = pwobj.Integer(kwargs.get('micId', None))
        self.micName = pwobj.String(kwargs.get('micName', None))

        # list particle x coordinate
        def _floatList(key):
            fl = pwobj.CsvList(pType=float)
            fl.set(kwargs.get(key, []))
            return fl

        self.x = _floatList('xCoord')
        # list particle y coordinate
        self.y = _floatList('yCoord')
        # list particle defocus
        self.defocus = _floatList('defocus')
        # list particle defocus difference
        self.defocusDiff = _floatList('defocusDiff')

        self.stdev = pwobj.Float()  # defocus stdev
        self.n = pwobj.Integer()  # number particles in the micrograph

    def addEntry(self, x, y, defocus, defocusDiff):
        self.x.append(x)
        self.y.append(y)
        self.defocus.append(defocus)
        self.defocusDiff.append(defocusDiff)

    def computeStats(self):
        self.n.set(len(self.defocus))
        self.stdev.set(np.std(self.defocus))

    def getCoordinates(self):
        return self.x, self.y

