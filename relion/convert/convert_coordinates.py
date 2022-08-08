# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] SciLifeLab, Stockholm University
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
# **************************************************************************
"""
Utility functions for conversions that will be used from both
newer Relion3.1 routines and old ones.
"""

import os
from emtable import Table
import logging
logger = logging.getLogger(__name__)

import pyworkflow.utils as pwutils
from pwem.constants import NO_INDEX


def openStar(fn, extraLabels=False):
    # We are going to write metadata directly to file to do it faster
    f = open(fn, 'w')
    s = """
data_

loop_
_rlnCoordinateX
_rlnCoordinateY
"""
    if extraLabels:
        s += "_rlnClassNumber\n"
        s += "_rlnAutopickFigureOfMerit\n"
        s += "_rlnAnglePsi\n"
    f.write(s)
    return f


def writeSetOfCoordinates(posDir, coordSet, getStarFileFunc, scale=1):
    """ Convert a SetOfCoordinates to Relion star files.
    Params:
        posDir: the output directory where to generate the files.
        coordSet: the input SetOfCoordinates that will be converted.
         getStarFileFunc: function object that receives the micrograph name
            and return the coordinates star file (only the base filename).
        scale: pass a value if the coordinates have a different scale.
            (for example when extracting from micrographs with a different
            pixel size than during picking)
    """

    # Create a dictionary with the pos filenames for each micrograph
    posDict = {}
    for mic in coordSet.iterMicrographs():
        starFile = getStarFileFunc(mic)
        if starFile is not None:
            posFn = os.path.basename(starFile)
            posDict[mic.getObjId()] = os.path.join(posDir, posFn)

    f = None
    lastMicId = None

    extraLabels = coordSet.getFirstItem().hasAttribute('_rlnClassNumber')
    doScale = abs(scale - 1) > 0.001

    for coord in coordSet.iterItems(orderBy='_micId'):
        micId = coord.getMicId()

        if micId != lastMicId:
            if micId not in posDict:
                logger.warning(f"Warning: micId {micId} not found")
                continue
            # we need to close previous opened file
            if f:
                f.close()
            f = openStar(posDict[micId], extraLabels)
            lastMicId = micId

        if doScale:
            x = coord.getX() * scale
            y = coord.getY() * scale
        else:
            x = coord.getX()
            y = coord.getY()

        if not extraLabels:
            f.write("%d %d \n" % (x, y))
        else:
            f.write("%d %d %d %0.6f %0.6f\n"
                    % (x, y,
                       coord._rlnClassNumber,
                       coord._rlnAutopickFigureOfMerit,
                       coord._rlnAnglePsi))

    if f:
        f.close()

    return posDict.values()


def writeSetOfCoordinatesXmipp(posDir, coordSet, ismanual=True, scale=1):
    """ Write a pos file on metadata format for each micrograph
    on the coordSet.
    Params:
        posDir: the directory where the .pos files will be written.
        coordSet: the SetOfCoordinates that will be read."""

    boxSize = coordSet.getBoxSize() or 100
    state = 'Manual' if ismanual else 'Supervised'

    # Create a dictionary with the pos filenames for each micrograph
    posDict = {}
    for mic in coordSet.iterMicrographs():
        micIndex, micFileName = mic.getLocation()
        micName = os.path.basename(micFileName)

        if micIndex != NO_INDEX:
            micName = '%06d_at_%s' % (micIndex, micName)

        posFn = os.path.join(posDir, pwutils.replaceBaseExt(micName, "pos"))
        posDict[mic.getObjId()] = posFn

    f = None
    lastMicId = None

    for coord in coordSet.iterItems(orderBy='_micId'):
        micId = coord.getMicId()

        if micId != lastMicId:
            # we need to close previous opened file
            if f:
                f.close()
            f = openMd(posDict[micId], state)
            lastMicId = micId
        if scale != 1:
            x = coord.getX() * scale
            y = coord.getY() * scale
        else:
            x = coord.getX()
            y = coord.getY()
        f.write(" %06d   1   %d  %d  %d   %06d\n"
                % (coord.getObjId(), x, y, 1, micId))

    if f:
        f.close()

    # Write config.xmd metadata
    configFn = os.path.join(posDir, 'config.xmd')
    writeCoordsConfig(configFn, boxSize, state)

    return posDict.values()


def writeCoordsConfig(configFn, boxSize, state):
    """ Write the config.xmd file needed for Xmipp picker.
    Params:
        configFn: The filename were to store the configuration.
        boxSize: the box size in pixels for extraction.
        state: picker state
    """
    # Write config.xmd metadata
    logger.debug(f"writeCoordsConfig: state={state}")
    table = Table(columns=['particleSize', 'pickingState'])
    table.addRow(int(boxSize), state)
    table.write(configFn, tableName='properties')


def openMd(fn, state='Manual'):
    # We are going to write metadata directly to file to do it faster
    f = open(fn, 'w')
    ismanual = state == 'Manual'
    block = 'data_particles' if ismanual else 'data_particles_auto'
    s = """# XMIPP_STAR_1 *
#
data_header
loop_
 _pickingMicrographState
%s
%s
loop_
 _itemId
 _enabled
 _xcoor
 _ycoor
 _cost
 _micrographId
""" % (state, block)
    f.write(s)
    return f


def writeMicCoordinates(mic, coordList, outputFn, getPosFunc=None):
    """ Write the pos file as expected by Xmipp with the coordinates
    of a given micrograph.
    Params:
        mic: input micrograph.
        coordList: list of (x, y) pairs of the mic coordinates.
        outputFn: output filename for the pos file .
        isManual: if the coordinates are 'Manual' or 'Supervised'
        getPosFunc: a function to get the positions from the coordinate,
            it can be useful for scaling the coordinates if needed.
    """
    if getPosFunc is None:
        getPosFunc = lambda coord: coord.getPosition()

    extraLabels = coordList[0].hasAttribute('_rlnAutopickFigureOfMerit')
    f = openStar(outputFn, extraLabels)

    for coord in coordList:
        x, y = getPosFunc(coord)
        if not extraLabels:
            f.write("%d %d \n" % (x, y))
        else:
            f.write("%d %d %d %0.6f %0.6f\n"
                    % (x, y,
                       coord._rlnClassNumber,
                       coord._rlnAutopickFigureOfMerit,
                       coord._rlnAnglePsi))

    f.close()
