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

import sys
import matplotlib as mpl
import numpy as np

from pwem.viewers.plotter import plt
from pwem.emlib.image import ImageHandler
from pyworkflow.viewer import ProtocolViewer

from .viewer_base import *
from ..objects import CtfRefineGlobalInfo
from ..protocols import ProtRelionCtfRefinement


class ProtCtfRefineViewer(ProtocolViewer):
    """ Viewer for Relion CTF refine results. """
    _label = 'ctf refine viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtRelionCtfRefinement]

    def __init__(self,  **kwargs):
        ProtocolViewer.__init__(self,  **kwargs)
        self.protocol._initialize()
        self._micInfoList = None
        self.xMax = None
        self.yMax = None
        self._currentMicIndex = 0
        self._oldCurrentMicIndex = 0
        self._currentMicId = 1
        self._loadAnalyzeInfo()

    def _defineParams(self, form):
        self._env = os.environ.copy()
        showBeamTilt = self.protocol.doBeamtiltEstimation.get()
        showTrefoil = self.protocol.doEstimateTrefoil.get() and showBeamTilt
        showTetrafoil = self.protocol.doEstimate4thOrder.get()
        showDefocus = self.protocol.doCtfFitting.get()
        showAnisoMag = self.protocol.estimateAnisoMag.get()

        form.addSection(label="Results")
        form.addParam('useMatplotlib', params.BooleanParam, default=True,
                      label='Use matplotlib for display',
                      condition="not {}".format(showDefocus),
                      help='If False, images will be displayed with ImageJ')

        form.addParam('displayAnisoMag', params.LabelParam,
                      label="Show X/Y mag. anisotropy estimation",
                      condition="{}".format(showAnisoMag),
                      help="Display four images (X and Y): (1) phase "
                           "differences from which this estimate was "
                           "derived and\n(2) the model fitted through it.")
        form.addParam('displayDefocus', params.LabelParam,
                      label="Show defocus estimation",
                      condition="{}".format(showDefocus),
                      help="Display the defocus estimation.\n "
                           "Plot defocus difference (Angstroms) vs "
                           "position in micrograph\n"
                           "You may move between micrographs by using: \n\n"
                           "Left/Right keys (move +1/-1 micrograph)\n"
                           "Up/Down keys (move +10/-10 micrographs)\n"
                           "Page_up/Page_down keys (move +100/-100 "
                           "micrographs)\n"
                           "Home/End keys (move +1000/-1000 micrographs)")
        form.addParam('displayBeamTilt', params.LabelParam,
                      label="Show beam tilt estimation",
                      condition="{} and not {}".format(showBeamTilt, showTrefoil),
                      help="Display two images: (1) phase differences from "
                           "which this estimate was derived and\n"
                           "(2) the model fitted through it.")
        form.addParam('displayTrefoil', params.LabelParam,
                      label="Show beam tilt and 3-fold astigmatism estimation",
                      condition="{}".format(showTrefoil),
                      help="Display two images: (1) phase differences from "
                           "which this estimate was derived and\n"
                           "(2) the model fitted through it.")
        form.addParam('displayTetrafoil', params.LabelParam,
                      label="Show 4-fold aberrations estimation",
                      condition="{}".format(showTetrafoil),
                      help="Display two images: (1) phase differences from "
                           "which this estimate was derived and\n"
                           "(2) the model fitted through it.")
        form.addParam('displayParticles', params.LabelParam,
                      label="Display particles",
                      help="See the particles with the new CTF "
                           "values")

    def _getVisualizeDict(self):
        return {
            'displayAnisoMag': self._displayAnisoMag,
            'displayDefocus': self._displayDefocus,
            'displayBeamTilt': self._displayBeamTilt,
            'displayTrefoil': self._displayTrefoil,
            'displayTetrafoil': self._displayTetrafoil,
            'displayParticles': self._displayParticles
        }

    def _displayAnisoMag(self, param=None):
        obs_x = self.protocol._getFileName("mag_obs_x", og=1)
        obs_y = self.protocol._getFileName("mag_obs_y", og=1)
        fit_x = self.protocol._getFileName("mag_fit_x", og=1)
        fit_y = self.protocol._getFileName("mag_fit_y", og=1)

        return self._showImages(obs_x, obs_y, fit_x, fit_y,
                                title="Anisotropy magnification")

    def _displayDefocus(self, e=None):
        """Show matplotlib with defocus values."""
        micInfo = self._micInfoList[self._currentMicIndex]
        self._currentMicId = micInfo.micId.get()
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

        self.plotter = EmPlotter(windowTitle="CTF Refinement", x=1, y=2)
        self._plotDefocusStdev()

        self.fig = self.plotter.getFigure()
        self.ax2 = self.plotter.createSubPlot(
            self._getTitle(micInfo), "Mic-Xdim", "Mic-Ydim",
            xpos=1, ypos=2)

        # call self.press after pressing any key
        self.fig.canvas.mpl_connect('key_press_event', self.press)

        # Maximize plot, valid for the 3 most common backends
        # Not sure if we need it since scipion installs it own TK
        # but I guess somebody may use the system one
        backend = mpl.get_backend()
        manager = plt.get_current_fig_manager()
        if backend == 'QT':
            # Option 1
            # QT backend
            manager.window.showMaximized()
        elif backend == 'TkAgg':
            # Option 2
            # TkAgg backend
            manager.resize(*manager.window.maxsize())
        elif backend == 'WX':
            # Option 3
            # WX backend
            manager.frame.Maximize(True)

        self.show()

    def _showImages(self, *imgs, **kwargs):
        title = kwargs.get('title', '')
        if self.useMatplotlib:
            return self._showImagesMatplotlib(title, *imgs)
        else:
            return [DataView(img) for img in imgs]

    def _showImagesMatplotlib(self, title, *imgs):
        ih = ImageHandler()
        xdim = 2 if len(imgs) > 2 else 1
        ydim = 2
        plotter = EmPlotter(windowTitle=title, x=xdim, y=ydim, figsize=(8, 6))
        positions = [(1, 1), (1, 2), (2, 1), (2, 2)]
        for i, imgFn in enumerate(imgs):
            x, y = positions[i]
            ax = plotter.createSubPlot("", "x", "y", x, y)
            img = ih.read(imgFn)
            ax.imshow(img.getData(), cmap='jet')

        return [plotter]

    def _displayBeamTilt(self, param=None):
        beamtilt_obs = self.protocol._getFileName("beamtilt_obs", og=1)
        beamtilt_fit = self.protocol._getFileName("beamtilt_fit", og=1)

        return self._showImages(beamtilt_obs, beamtilt_fit,
                                title="Beam Tilt")

    def _displayTrefoil(self, param=None):
        trefoil_obs = self.protocol._getFileName("beamtilt_obs", og=1)
        trefoil_fit = self.protocol._getFileName("trefoil_fit", og=1)

        return self._showImages(trefoil_obs, trefoil_fit,
                                title="Trefoil")

    def _displayTetrafoil(self, param=None):
        tetrafoil_obs = self.protocol._getFileName("tetrafoil_obs", og=1)
        tetrafoil_fit = self.protocol._getFileName("tetrafoil_fit", og=1)

        return self._showImages(tetrafoil_obs, tetrafoil_fit,
                                title="Tetrafoil")

    def createScipionPartView(self, filename):
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

    def _displayParticles(self, param=None):
        views = []
        fn = self.protocol.outputParticles.getFileName()
        v = self.createScipionPartView(fn)
        views.append(v)
        return views

    # ------------------- UTILS functions -------------------------
    def _loadAnalyzeInfo(self):
        # Only load once
        if self._micInfoList is None:
            ctfInfoFn = self.protocol._getFileName("ctf_sqlite")
            if not os.path.exists(ctfInfoFn):
                ctfInfo = self.protocol.createGlobalInfo(ctfInfoFn)
            else:
                ctfInfo = CtfRefineGlobalInfo(ctfInfoFn)
            self._micInfoList = [mi.clone() for mi in ctfInfo]
            self.xMax, self.yMax = ctfInfo.getMaxXY()
            self.ctfInfoMapper = ctfInfo
            ctfInfo.close()
            self.len_micInfoList = len(self._micInfoList)
            micList = [mi.micId.get() for mi in self._micInfoList]
            # instead of creating this dict we have
            # access the mapper and make a query to the database
            self.micDict = {k: v for v, k in enumerate(micList)}

    def onClick(self, event):
        # try is needed because if clicked outside plot
        # xdata, ydata are Nonetype
        try:
            ix, iy = int(round(event.xdata)), int(round(event.ydata))
            # once the user has selected a point
            # he have a pair of float numbers,
            # search for the closest micrograph
            # with the right stdev in a neighbourhood
            if ix <= self.maxMicId:
                while ix not in self.micDict:
                    ix += 1
                iix = self.micDict[ix]
                start = max(0, iix-40)
                end = min(iix+40, self.maxMicId)
                dist = np.sqrt((np.array(self.x[start:end]) - ix) ** 2 +
                               (np.array(self.y[start:end]) - iy) ** 2)
                self._currentMicId = self.x[start + np.argmin(dist)]
                self._oldCurrentMicIndex = self._currentMicIndex
                self._currentMicIndex = self.micDict[self._currentMicId]
                self.show()
        except:
            pass

    def _plotDefocusStdev(self, e=None):
        self.fig = self.plotter.getFigure()
        self.ax1 = self.plotter.createSubPlot("Defocus stdev per Micrograph\n"
                                              "Click on any point to get the "
                                              "corresponding micrograph\n "
                                              "in the defocus plot",
                                              "# Micrograph", "stdev",
                                              xpos=1, ypos=1)
        self.fig.canvas.mpl_connect('button_press_event', self.onClick)
        self.ax1.grid(True)
        self.maxMicId = self._micInfoList[-1].micId.get()
        self.x = [mi.micId.get() for mi in self._micInfoList]
        self.y = [mi.stdev.get() for mi in self._micInfoList]
        self.ax1.scatter(self.x, self.y, s=50, marker='o',
                         c='blue')

    def _getTitle(self, micInfo):
        return ("Use arrows or Page up/Down or Home/End to navigate.\n"
                "Mic = %s (%d)\nColorBar indicates defocus difference" %
                (micInfo.micName.get()[-40:], micInfo.micId))

    def press(self, event):
        """ Change the currently shown micrograph
        when a key is pressed (increment/decrement)
        """
        sys.stdout.flush()

        if event.key == 'q':
            plt.close('all')

        shiftDict = {
            'left': -1, 'right': 1,
            'up': 10, 'down': -10,
            'pageup': -100, 'pagedown': 100,
            'home': -1000, 'end': 1000
        }
        # if pressed key is not left, up, etc, do nothing
        if event.key in shiftDict:
            shift = shiftDict[event.key]
            # Check the new micrograph index is between first and last
            newIndex = self._currentMicIndex + shift
            self._oldCurrentMicIndex = self._currentMicIndex
            self._currentMicIndex = min(max(0, newIndex),
                                        self.len_micInfoList - 1)
            self.show()

    def show(self, event=None):
        """ Draw plot """
        # stdev plot
        self.ax1.scatter(self.x[self._oldCurrentMicIndex],
                         self.y[self._oldCurrentMicIndex], s=50,
                         marker='o', c='blue')
        self.ax1.scatter(self.x[self._currentMicIndex],
                         self.y[self._currentMicIndex], s=50,
                         marker='o', c='red')

        # defocus plot
        micInfo = self._micInfoList[self._currentMicIndex]
        if event is None:
            # I need to clear the plot otherwise
            # old points are not removed
            self.ax2.clear()
            self.ax2.margins(0.05)
            self.ax2.set_title(self._getTitle(micInfo))
            newFontSize = self.plotter.plot_axis_fontsize + 2
            self.ax2.set_xlabel("Mic Xdim (px)", fontsize=newFontSize)
            self.ax2.set_ylabel("Mic Ydim (px)", fontsize=newFontSize)

            self.ax2.set_xlim(0, self.xMax)  # np.max(micInfo.x))
            self.ax2.set_ylim(0, self.yMax)  # np.max(micInfo.y))
            self.ax2.grid(True)

            # if I do not use subplots_adjust the window shrinks
            #  after redraw
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        sc2 = self.ax2.scatter(micInfo.x, micInfo.y,
                               c=micInfo.defocusDiff, s=100, marker='o')
        self.plotter.getColorBar(sc2)
        self.plotter.show()
