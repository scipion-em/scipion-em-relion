import sqlite3
from pyworkflow.protocol.params import LabelParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pyworkflow.em.data import SetOfParticles
from relion.protocols.protocol_ctf_refinement import ProtRelionCtfRefinement
from itertools import izip
from pyworkflow.em.viewers.plotter import plt, EmPlotter
import matplotlib as mpl
import pyworkflow.em.viewers.showj as showj
from pyworkflow.em.viewers import ObjectView, DataView
import sys
import os


class ProtCtfREfineViewer(ProtocolViewer):
    """ viewer of relaion cTF refine"""
    _label = 'CTF Refine Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtRelionCtfRefinement]

    def __init__(self,  **kwargs):
        ProtocolViewer.__init__(self,  **kwargs)
        # TODO: move this inizialization to protocol
        # create temporary database to store values
        self.step = 1  # next micrography

        def createDB(tableName):
            """ create database and table"""
            # TODO: IS this sectio is moved to protocol
            # data base must be written in disk rather than in memory
            conn = sqlite3.connect(":memory:")
            # conn = sqlite3.connect(self._getExtraPath("kk.sqlite3"))
            c = conn.cursor()
            commandCreateTable = """
                    CREATE TABLE {tableName}(
                         id           integer primary key autoincrement,
                         micId        int,
                         micName      text,
                         coordX       float,
                         coordY       float,
                         defocusDiff  float
                         );"""
            c.execute(commandCreateTable.format(tableName=tableName))
            return c, conn

        self.tableName = 'defoci'

        self.c, self.conn = createDB(self.tableName)
        self.doDb = True

    def _defineParams(self, form):
        self._env = os.environ.copy()
        showBeamTilt = self.protocol.doBeamtiltEstimation.get()
        form.addSection(label="Defocus")
        form.addParam('displayDefocus', LabelParam,
                      label="Show Defocus",
                      help="Display the defocus estimation.\n "
                           "Plot defocus difference vs positin in micrograph\n"
                           "You may move between microgaphs by typing:\:"
                           "left/right keys (move +1/-1 micrograph)\n"
                           "up/down keys (move +10/-10 micrographs)\n"
                           "page_up/page_down keys (move +100/-100 "
                           "micrographs)\n"
                           "home/end keys (move +1000/-1000 micrographs)\n")
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
            'displayParticles': self._displayParticles
        }

    def _visualizeDefocus(self, e=None):
        """Show matplotlib with defocus values."""
        # read input and output metadata
        if self.doDb:  # TODO MOVE DATABASE INSERTIONS TO  PROTOCOL
            outPartSet = SetOfParticles(
                filename=self.protocol._getPath("particles.sqlite"))
            inPartSet = SetOfParticles(
                filename=self.protocol.inputParticles.get().getFileName())
            # compute difference in defocus
            # and save in database
            data = []
            resultDict = {}
            resultDict['tableName'] = self.tableName

            for p1, p2 in izip(inPartSet, outPartSet):
                p1D = (p1.getCTF().getDefocusU() +
                       p1.getCTF().getDefocusV())/2.0
                p2D = (p2.getCTF().getDefocusU() +
                       p2.getCTF().getDefocusV())/2.0
                coordinate = p1.getCoordinate()
                micId = coordinate.getMicId()
                micName = coordinate.getMicName()
                coordX = coordinate.getX()
                coordY = coordinate.getY()
                defocusDiff = p2D - p1D
                data.append((micId, micName, coordX, coordY,  defocusDiff))
            try:
                sql = "INSERT INTO %s  " \
                      " (micID, micName, coordX, coordY, defocusDiff)" \
                      " VALUES " \
                      " (?, ?, ?, ?, ?)" % self.tableName
                self.c.executemany(sql, data)
            except sqlite3.IntegrityError as e:
                print('sqlite error: ', e.args[0])  # something went wrong

            self.conn.commit()

            # create index on micId
            sql = "CREATE INDEX micIDindex ON %s(micId);" % self.tableName
            self.c.execute(sql)

            # find first and last micId
            sql = "SELECT min(micId), max(micId) FROM %s" % self.tableName
            self.c.execute(sql)
            row = self.c.fetchone()
            self.smallerMicId = int(row[0])
            self.higherMicId = int(row[1])
            self.showMicWitID = self.smallerMicId

            # get micrograph name
            sql = "SELECT micName FROM %s WHERE micID=%s LIMIT 1" %\
                  (self.tableName, self.smallerMicId)
            self.c.execute(sql)
            row = self.c.fetchone()
            self.micName = row[0]

            self.doDb = False

            # find plot dimensions
            # I would rather use the micrograph dimensions but I do
            # not know how to get them
            sql = "SELECT max(coordX), max(coordY) FROM %s" % self.tableName
            self.c.execute(sql)
            row = self.c.fetchone()
            self.xSize = int(row[0])
            self.ySize = int(row[1])

            # disable default binding for arrows
            # because I want to use them
            # to navigate between micrographs
            mpl.rcParams['keymap.back'].remove('left')
            mpl.rcParams['keymap.forward'].remove('right')

        self.plotter = EmPlotter(windowTitle="CTF Refinement")

        self.fig = self.plotter.getFigure()
        self.ax = \
            self.plotter.createSubPlot(self._getTitle(),
                                       "Mic-Xdim",
                                       "Mic-Ydim",
                                       projection='3d')
        # call self.press after pressing any key
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        # scatter plots loss colormap after rotation, so
        # if the user rotates the canvas repaint.
        # I guess this is a bug in matplotlib
        self.fig.canvas.mpl_connect('draw_event', self.show)

        self.show()

    def _getTitle(self):
        return ("use arrows or page up/down or home/end to navigate.\n"
                "Displaying Mic = %s (%d)" %
                (self.micName, self.showMicWitID))

    def getData(self):
        sql = "SELECT coordX, coordY, defocusDiff, micName " \
              "FROM %s " \
              "WHERE micId = %d"
        while True:
            sqlComamnd = sql % (self.tableName,  self.showMicWitID)
            self.c.execute(sqlComamnd)
            rows = self.c.fetchall()
            if len(rows) == 0 \
                    and self.showMicWitID >= self.smallerMicId \
                    and self.showMicWitID <= self.higherMicId:
                self.showMicWitID += self.step
            else:
                break


        x = [item[0] for item in rows]
        y = [item[1] for item in rows]
        defocus = [item[2] for item in rows]
        micName = rows[0][3]
        return x, y, defocus, micName

    def press(self, event):
        """ if a key is pressed increment/decrement
        the showMicWitID"""
        sys.stdout.flush()
        if event.key == "left":
            self.showMicWitID -= 1
            self.step = -1
        elif event.key == "right":
            self.showMicWitID += 1
            self.step = +1
        elif event.key == "up":
            self.showMicWitID -= 10
            self.step = -1
        elif event.key == "down":
            self.showMicWitID += 10
            self.step = +1
        elif event.key == "pageup":
            self.showMicWitID -= 100
            self.step = -1
        elif event.key == "pagedown":
            self.showMicWitID += 100
            self.step = +1
        elif event.key == "home":
            self.showMicWitID -= 1000
            self.step = -1
        elif event.key == "end":
            self.showMicWitID += 1000
            self.step = +1
        elif event.key == "q":
            quit()
        # else:
        #    print('press', event.key)

        # you should no go beyond the last micrograph
        # and before the first one
        if self.showMicWitID < self.smallerMicId:
            self.showMicWitID = self.smallerMicId
        elif self.showMicWitID > self.higherMicId:
            self.showMicWitID = self.higherMicId

        self.show()

    def show(self, event=None):
        """ Draw plot """
        if event is None:
            self.x, self.y, self.defocus, self.micName = \
                self.getData()
            self.ax.clear()
            self.ax.margins(0.05)
            self.ax.set_title(self._getTitle())
            self.ax.set_xlabel("Mic Xdim (px)",
                               fontsize=self.plotter.plot_axis_fontsize+2)
            self.ax.set_ylabel("Mic Ydim (px)",
                               fontsize=self.plotter.plot_axis_fontsize+2)
            self.ax.set_zlabel("Defocus Difference (A)",
                               fontsize=self.plotter.plot_axis_fontsize+2)
            self.ax.grid(True)

            self.ax.set_xlim(0, self.xSize)
            self.ax.set_ylim(0, self.ySize)

            # if I do not use subplots_adjust the window shrinks
            #  after redraw
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

#        im = self.plotter.plotScatter(self.x, self.y, zs=self.defocus,
        # IT will be niver to create a plotScatter function in scipion but
        # scipion is frozen
        im = self.ax.scatter(self.x, self.y, zs=self.defocus,
                             c=self.defocus,
                             s=100,
                             marker='o')

        self.plotter.getColorBar(im)
        self.plotter.show()

    def _displayBeamTilt(self, paramName=None):
        phaseDifferenceFn = self.protocol.fileWithPhaseDifferenceName()
        modelFitFn = self.protocol.fileWithModelFitterName()
        # TODO: how can I change the window title?

        return [DataView(phaseDifferenceFn), DataView(modelFitFn)]

    def createScipionPartView(self, filename, viewParams={}):
        inputParticlesId = self.protocol.inputParticles.get().strId()
        labels = 'enabled id _size _filename ' \
                 ' _ctfModel._defocusU _ctfModel._defocusV '
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
