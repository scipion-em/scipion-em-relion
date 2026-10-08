# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow.object as pwobj
from pyworkflow.protocol.constants import STATUS_NEW
import pyworkflow.utils as pwutils
from pwem.protocols import ProtParticlePickingAuto

import relion.convert as convert
from .protocol_base import ProtRelionBase
from .protocol_streaming_base import RelionStreamingBase

# Module level: these helpers are called unbound on light test harnesses.
PICKING_STEP_NAMES = ('pickMicrographStep', 'pickMicrographListStep')


class ProtRelionAutopickBase(RelionStreamingBase, ProtParticlePickingAuto,
                             ProtRelionBase):
    """ Base class for auto-picking protocols in Relion.
    """
    _label = None

    def _loadSet(self, inputSet, SetClass, getKeyFunc, watermarkAttr=None,
                 knownIds=None):
        """Discover the items added to one input stream since the last poll.

        pwem rebuilds the Set from its storage filename and walks it whole.
        This queries the ids above that stream's watermark and loads only
        those, so a poll costs what just arrived rather than everything
        the stream has produced so far.
        """
        watermarkAttr = watermarkAttr or '_lastInputId'

        if knownIds is None:
            knownIds = self._getKnownStreamIds(watermarkAttr)

        newItems, producerClosed, terminalConsistent = (
            self._discoverNewInputItems(inputSet, watermarkAttr, knownIds))

        newItemDict = {}

        for item in newItems:
            itemId = item.getObjId()

            if itemId in knownIds:
                continue

            knownIds.add(itemId)
            itemKey = getKeyFunc(item)

            if itemKey not in self.micDict:
                newItemDict[itemKey] = item

        return newItemDict, producerClosed and terminalConsistent

    def _loadMics(self, micSet):
        """Micrographs discovered on their own watermark."""
        return self._loadSet(micSet, None, lambda mic: mic.getMicName(),
                             watermarkAttr='_lastMicId')

    def _loadCTFs(self, ctfSet):
        """CTFs discovered on a separate watermark: they are their own
        stream and advance independently of the micrographs."""
        return self._loadSet(ctfSet, None,
                             lambda ctf: ctf.getMicrograph().getMicName(),
                             watermarkAttr='_lastCtfId')

    def _checkNewInput(self):
        # Refresh logical input state directly. Do not gate discovery on
        # storage filenames or filesystem modification times.
        micDict, self.streamClosed = self._loadInputList()
        outputStep = self._getFirstJoinStep()

        if micDict:
            deps = self._insertNewMicsSteps(micDict.values())
            if outputStep is not None:
                outputStep.addPrerequisites(*deps)
            self.updateSteps()

    def _getFinishedPickingMicNames(self):
        """Micrograph names carried by finished picking steps."""
        return self._collectStepArgKeys(PICKING_STEP_NAMES, keyType=str)

    def _getPublishedPickingMicIds(self):
        """Micrograph ids already represented in the output coordinates.

        An id query on the output, not a walk over it.
        """
        micIds = self._getOutputUniqueValues(
            getattr(self, 'outputCoordinates', None), '_micId')

        return set() if micIds is None else micIds

    def _checkNewOutput(self):
        """Publish finished picking without DONE sidecars.

        pwem decides this with extra/DONE/mic_*.TXT plus DONE/all.TXT.
        The persisted step graph already records what finished and the
        output Set what was published, so neither file is consulted and a
        cleaned extra/ cannot make finished work look pending.
        """
        if getattr(self, 'finished', False):
            return

        publishedMicIds = self._getPublishedPickingMicIds()
        finishedMicNames = self._getFinishedPickingMicNames()

        listOfMics = list(self.micDict.values())
        newDone = [mic for mic in listOfMics
                   if mic.getMicName() in finishedMicNames
                   and mic.getObjId() not in publishedMicIds]

        doneCount = len([mic for mic in listOfMics
                         if mic.getObjId() in publishedMicIds])
        allDone = doneCount + len(newDone)

        self.debug('_checkNewOutput: mics=%s published=%s newDone=%s'
                   % (len(listOfMics), doneCount, len(newDone)))

        self.finished = self.streamClosed and allDone == len(listOfMics)
        streamMode = (pwobj.Set.STREAM_CLOSED
                      if self.finished else pwobj.Set.STREAM_OPEN)

        if newDone:
            self._updateOutputCoordSet(newDone, streamMode)
        elif not self.finished:
            self._streamingSleepOnWait()
            return

        if self.finished:
            self._updateStreamState(streamMode)

            # The lifecycle here is still the classic one: createOutputStep
            # was scheduled with wait=True and only runs once something
            # releases it. Dropping this leaves it WAITING for good and the
            # protocol never finishes, however complete the output is.
            outputStep = self._getFirstJoinStep()

            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

    def _pickMicrograph(self, mic, *args):
        """ This method should be invoked only when working in streaming mode.
        """
        self._pickMicrographList([mic], *args)

    def _pickMicrographList(self, micList, *args):
        if not micList:
            return

        micsDir = self._createTmpMicsDir(micList)
        micStar = os.path.join(micsDir, 'input_micrographs.star')
        writer = convert.createWriter(rootDir=micsDir, outputDir=micsDir)
        writer.writeSetOfMicrographs(micList, micStar)
        try:
            # pickMicrographListStep (pwem) has no exception boundary of
            # its own around this call - a single relion_autopick crash
            # for the whole batch would otherwise propagate uncaught and
            # abort the entire streaming run instead of being reported
            # like every other missing-output case below.
            self._pickMicrographsFromStar(micStar, micsDir, *args)
        except Exception as e:
            self.error(
                "ERROR: Autopick failed for micrograph batch starting at "
                "%s with the exception %s" % (micList[0].getObjId(), e)
            )

        # Move each expected coordinates file explicitly. Do not let a
        # successful wildcard move hide a missing output for one mic in a batch.
        for mic in micList:
            fileName = 'mic_%06d_autopick.star' % mic.getObjId()
            srcFile = os.path.join(micsDir, fileName)
            dstFile = self._getTmpPath(fileName)

            if os.path.exists(srcFile):
                pwutils.moveFile(srcFile, dstFile)
            elif not os.path.exists(dstFile):
                self.warning(
                    "Missing autopick output for micrograph %s: %s / %s"
                    % (mic.getObjId(), srcFile, dstFile)
                )

    def _createSetOfCoordinates(self, micSet, suffix=''):
        """ Override this method to set the box size. """
        coordSet = ProtParticlePickingAuto._createSetOfCoordinates(
            self, micSet, suffix=suffix)
        coordSet.setBoxSize(self.getBoxSize())
        return coordSet

    def readCoordsFromMics(self, workingDir, micList, coordSet):
        """ Parse back the output star files and populate the SetOfCoordinates.
        """
        template = self._getTmpPath("mic_%06d_autopick.star")
        starFiles = [template % mic.getObjId() for mic in micList]
        convert.readSetOfCoordinates(coordSet, starFiles, micList)

    def _pickMicrographsFromStar(self, micStar, micsDir, *args):
        """ Should be defined in subclasses. """
        pass

    def getBoxSize(self):
        """ Return a reasonable box-size in pixels. """
        return None

    def getInputMicrographsPointer(self):
        return self.inputMicrographs

    def getInputMicrographs(self):
        return self.getInputMicrographsPointer().get()

    def __getMicListPrefix(self, micList):
        n = len(micList)
        if n == 0:
            raise ValueError("Empty micrographs list!")
        micsPrefix = 'mic_%06d' % micList[0].getObjId()
        if n > 1:
            micsPrefix += "-%06d" % micList[-1].getObjId()
        return micsPrefix

    def _createTmpMicsDir(self, micList):
        """ Create a temporary path to work with a list of micrographs. """
        micsDir = self._getTmpPath(self.__getMicListPrefix(micList))
        pwutils.makePath(micsDir)
        return micsDir
