# **************************************************************************
# *
# * Relion particle extraction in streaming loses no micrograph and no
# * particle, never publishes a micrograph by halves, and always ends.
# *
# **************************************************************************

import os
import tempfile
import unittest
from collections import OrderedDict

from pwem.objects import Coordinate, Particle

from relion.protocols.protocol_extract_particles import (
    ProtRelionExtractParticles)



class _InMemoryFailedExtraction:
    """The failed-micrograph list, kept in memory instead of in extra/."""

    def _failedExtractMicIds(self):
        return set(self.__dict__.setdefault('failedIds', set()))

    def _recordFailedExtraction(self, mic):
        self.__dict__.setdefault('failedIds', set()).add(mic.getObjId())

    def _forgetFailedExtraction(self, micIds):
        self.__dict__.setdefault('failedIds', set()).difference_update(micIds)

class _Mic:
    def __init__(self, objId):
        self._objId = objId

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return "mic_%03d" % self._objId

    def getAttributeValue(self, name, default=None):
        return default

    def getCTF(self):
        return None


def _coord(micId, x, y, objId):
    coord = Coordinate()
    coord.setPosition(x, y)
    coord.setMicId(micId)
    coord.setObjId(objId)
    return coord


class _CoordSet:
    """Coordinates that turn up a poll at a time."""

    def __init__(self, coords=(), closed=False):
        self.coords = list(coords)
        self.closed = closed

    def loadAllProperties(self):
        pass

    def isStreamClosed(self):
        return self.closed

    def getUniqueValues(self, attribute, where=None):
        last = int(where.split('>')[1]) if where else 0
        values = [c.getObjId() if attribute == 'id' else c.getMicId()
                  for c in self.coords if c.getObjId() > last]
        return sorted(set(values))

    def iterItems(self, where=None):
        micId = int(where.split('=')[1])
        return [c for c in self.coords if c.getMicId() == micId]


class _Extract(_InMemoryFailedExtraction, ProtRelionExtractParticles):
    def __init__(self, coordSet=None, batchSize=1):
        self.micDict = OrderedDict()
        self.coordDict = {}
        self.initialIds = []
        self.streamClosed = False
        self._coordSet = coordSet
        self._batchSize = batchSize
        self.batches = []
        self.errors = []

    def getCoords(self):
        return self._coordSet

    def _getStreamingBatchSize(self):
        return self._batchSize

    def _getExtractArgs(self):
        return []

    def _insertExtractMicrographListStep(self, micList, prerequisites, *args):
        self.batches.append([mic.getObjId() for mic in micList])
        return len(self.batches)

    def _insertExtractMicrographStep(self, mic, prerequisites, *args):
        return self._insertExtractMicrographListStep([mic], prerequisites)

    def error(self, message):
        self.errors.append(message)

    def warning(self, message):
        pass

    def debug(self, message):
        pass


class TestNoMicrographIsLost(unittest.TestCase):
    def test_TheRemainderOfABatchWaitsForTheNextPoll(self):
        protocol = _Extract(batchSize=2)

        protocol._insertNewMicsSteps([_Mic(1), _Mic(2), _Mic(3)])
        protocol._insertNewMicsSteps([_Mic(4)])

        self.assertEqual([[1, 2], [3, 4]], protocol.batches)

    def test_AMicrographSeenBeforeItsCoordinatesWaitsForThem(self):
        coords = _CoordSet()
        protocol = _Extract(coords)

        self.assertEqual({}, protocol._loadInputCoords({"mic_001": _Mic(1)}))

        coords.coords.append(_coord(1, 10, 20, objId=1))
        ready = protocol._loadInputCoords({})

        self.assertEqual(["mic_001"], list(ready))
        self.assertEqual(1, len(protocol.coordDict[1]))

    def test_AMicrographWithNoCoordinatesAtAllIsLetGoWhenTheyClose(self):
        coords = _CoordSet()
        protocol = _Extract(coords)
        protocol._loadInputCoords({"mic_001": _Mic(1)})

        coords.closed = True

        self.assertEqual({}, protocol._loadInputCoords({}))
        self.assertFalse(protocol._micsWaitingForCoords())


class _Acq:
    def getMagnification(self):
        return 50000

    def getVoltage(self):
        return 300

    def getAmplitudeContrast(self):
        return 0.1

    def getSphericalAberration(self):
        return 2.7


class _Mics:
    def getAcquisition(self):
        return _Acq()


class _Output:
    def __init__(self):
        self.particles = []

    def append(self, particle):
        self.particles.append(particle.clone())


class _Reading(_Extract):
    def __init__(self, root, rows, scale=2.0):
        super().__init__()
        self._root = root
        self._scale = scale
        self._getPos = lambda c: (int(c.getX()), int(c.getY()))
        os.makedirs(os.path.join(root, 'tmp'))
        os.makedirs(os.path.join(root, 'extra'))
        with open(os.path.join(root, 'tmp', 'mic_000001_extract.star'), 'w') as fh:
            fh.write("data_\n\nloop_\n_rlnCoordinateX #1\n_rlnCoordinateY #2\n"
                     "_rlnImageName #3\n")
            fh.writelines("%s\n" % row for row in rows)
        open(os.path.join(root, 'tmp', 'mic_000001.mrcs'), 'w').close()

    def _getTmpPath(self, *paths):
        return os.path.join(self._root, 'tmp', *paths)

    def _getExtraPath(self, *paths):
        return os.path.join(self._root, 'extra', *paths)

    def getInputMicrographs(self):
        return _Mics()

    def getBoxScale(self):
        return self._scale


class TestParticlesAreReadExactlyOnce(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self._tmp.cleanup()

    def test_ReadingAgainAfterAFailedWriteDoesNotScaleTwice(self):
        protocol = _Reading(self._tmp.name, ["10 20 1@a.mrcs", "30 40 2@a.mrcs"])
        protocol.coordDict[1] = [_coord(1, 10, 20, 11), _coord(1, 30, 40, 12)]

        first, again = _Output(), _Output()
        protocol.readPartsFromMics([_Mic(1)], first)
        protocol.readPartsFromMics([_Mic(1)], again)

        self.assertEqual([(20, 40), (60, 80)],
                         [p.getCoordinate().getPosition() for p in first.particles])
        self.assertEqual([(20, 40), (60, 80)],
                         [p.getCoordinate().getPosition() for p in again.particles])

    def test_AMicrographThatFailsHalfWayPublishesNothing(self):
        protocol = _Reading(self._tmp.name, ["10 20 1@a.mrcs", "not-a-number 40 2@a.mrcs"])
        protocol.coordDict[1] = [_coord(1, 10, 20, 11), _coord(1, 30, 40, 12)]

        output = _Output()
        protocol.readPartsFromMics([_Mic(1)], output)

        self.assertEqual([], output.particles)
        self.assertTrue(protocol.errors)
        self.assertEqual({1}, protocol._failedExtractMicIds())


class _Closing(_Extract):
    def __init__(self, failed):
        super().__init__()
        self.micDict = OrderedDict([("mic_001", _Mic(1))])
        self.micsClosed = self.ctfsClosed = self.coordsClosed = True
        self._failed = failed
        self.published = []

    def _getPublishedExtractMicIds(self):
        return set()

    def _getFinishedExtractMicNames(self):
        return {"mic_001"}

    def _updateOutputPartSet(self, micList, streamMode):
        self.published.append([mic.getObjId() for mic in micList])
        if self._failed and micList:
            self._recordFailedExtraction(micList[0])

    def _getFirstJoinStep(self):
        return None

    def _streamingSleepOnWait(self):
        pass


class TestExtractionAlwaysEnds(unittest.TestCase):
    def test_AMicrographWithoutParticlesIsReadOnce(self):
        # Nothing of it reaches the output, so "not published yet" alone
        # would read it, and fail, again on every poll.
        protocol = _Closing(failed=True)
        protocol.micsClosed = False
        protocol.outputParticles = _Sized(5)

        ProtRelionExtractParticles._checkNewOutput(protocol)
        ProtRelionExtractParticles._checkNewOutput(protocol)

        self.assertEqual([[1]], protocol.published)

    def test_TheRunEndsWhenEveryMicrographWasHandled(self):
        protocol = _Closing(failed=False)
        protocol.outputParticles = _Sized(5)

        ProtRelionExtractParticles._checkNewOutput(protocol)

        self.assertTrue(protocol.finished)

    def test_EveryMicrographFailedFailsTheRun(self):
        protocol = _Closing(failed=True)

        with self.assertRaisesRegex(Exception, "Extraction failed for every"):
            ProtRelionExtractParticles._checkNewOutput(protocol)


class _Sized:
    def __init__(self, size):
        self._size = size

    def getSize(self):
        return self._size


if __name__ == '__main__':
    unittest.main()
