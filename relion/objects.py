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

import numpy as np

import pyworkflow.object as pwobj


# --------- Protocol CtfRefine Visualization Classes  -----------------

class CtfRefineGlobalInfo:
    """ Simple class to store visualization information related
    to Micrographs and Particles CTF information after
    Relion - ctf refinement protocol.
    """

    def __init__(self, filename):
        # The classes dict needs to be updated to register local objects
        classesDict = dict(pwobj.__dict__)
        classesDict['CtfRefineMicInfo'] = CtfRefineMicInfo
        self._infoSet = pwobj.Set(filename, classesDict=classesDict)

    def addMicInfo(self, micId, x, y, defocus, defocusDiff):
        pass

    def loadFromParticles(self, inputParts, outputParts):
        # compute difference in defocus and save in database
        micInfo = None
        lastMicId = None
        infoList = []
        # Coordinates of the particle with higher x and y values
        # These values are used for plotting
        xMax = 0
        yMax = 0

        def _avgDefocus(p):
            dU, dV, _ = p.getCTF().getDefocus()
            return (dU + dV) / 2.0

        for p1, p2 in zip(inputParts.iterItems(orderBy=['_micId', 'id']),
                          outputParts.iterItems(orderBy=['_micId', 'id'])):
            coord = p1.getCoordinate()
            micId = coord.getMicId()

            if micId != lastMicId:
                micInfo = CtfRefineMicInfo(micId=micId,
                                           micName=coord.getMicName())
                infoList.append(micInfo)
                lastMicId = micId

            p1D = _avgDefocus(p1)
            p2D = _avgDefocus(p2)
            x, y = coord.getPosition()
            if xMax < x:
                xMax = x
            if yMax < y:
                yMax = y
            micInfo.addEntry(x, y, p2D, p2D - p1D)

        for micInfo in infoList:
            micInfo.computeStats()
            self._infoSet.append(micInfo)

        self._infoSet._xMax = pwobj.Integer(xMax)
        self._infoSet._yMax = pwobj.Integer(yMax)
        self._infoSet.write()

    def getMaxXY(self):
        """ Return maximum value of coordinates"""
        mapper = self._infoSet._getMapper()
        # TODO: eventually remove this if
        # since in the future all instances will
        # have _xMax (mar 22nd 2019)
        if mapper.hasProperty('_xMax'):
            return int(mapper.getProperty('_xMax')), \
                   int(mapper.getProperty('_yMax'))
        else:
            return -1, -1

    def __iter__(self):
        for micInfo in self._infoSet:
            yield micInfo

    def close(self):
        self._infoSet.close()


class CtfRefineMicInfo(pwobj.Object):
    def __init__(self, **kwargs):
        pwobj.Object.__init__(self, **kwargs)

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
