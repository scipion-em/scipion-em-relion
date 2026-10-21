# **************************************************************************
# *
# * Relion auto-picking in streaming: a micrograph is published whole or
# * not at all, read once per run, and the run fails only when every
# * micrograph failed.
# *
# **************************************************************************

import os
import tempfile
import unittest
from collections import OrderedDict
from unittest import mock

from relion.protocols.protocol_autopick import ProtRelionAutopickBase


class _Mic:
    def __init__(self, objId):
        self._objId = objId

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return "mic_%03d" % self._objId


class _Coords:
    def __init__(self):
        self.items = []

    def append(self, coord):
        self.items.append(coord)


class _Picking(ProtRelionAutopickBase):
    def __init__(self, root):
        self._root = root
        self.errors = []
        self.micDict = OrderedDict()

    def _getTmpPath(self, *paths):
        return os.path.join(self._root, *paths)

    def _getExtraPath(self, *paths):
        os.makedirs(os.path.join(self._root, 'extra'), exist_ok=True)
        return os.path.join(self._root, 'extra', *paths)

    def error(self, message):
        self.errors.append(message)

    def warning(self, message):
        pass

    def debug(self, message):
        pass


def _readHalfOfTheSecond(mic, fileName, coordSet):
    coordSet.append((mic.getObjId(), 'first'))
    if mic.getObjId() == 2:
        raise ValueError("corrupted row")
    coordSet.append((mic.getObjId(), 'second'))


class TestAMicrographIsPublishedWhole(unittest.TestCase):
    def test_AStarThatFailsHalfWayPublishesNothingOfItsMicrograph(self):
        with tempfile.TemporaryDirectory() as root:
            protocol = _Picking(root)
            for micId in (1, 2):
                open(protocol._getTmpPath('mic_%06d_autopick.star' % micId), 'w').close()
            coords = _Coords()

            with mock.patch('relion.protocols.protocol_autopick.convert.readCoordinates',
                            _readHalfOfTheSecond):
                protocol.readCoordsFromMics(None, [_Mic(1), _Mic(2)], coords)

            self.assertEqual([(1, 'first'), (1, 'second')], coords.items)
            self.assertEqual({2}, protocol._failedPickingMicIds())

    def test_AFailedMicrographIsNotRead(self):
        with tempfile.TemporaryDirectory() as root:
            protocol = _Picking(root)
            protocol._recordFailedPicking(_Mic(1))
            coords = _Coords()

            with mock.patch('relion.protocols.protocol_autopick.convert.readCoordinates') as read:
                protocol.readCoordsFromMics(None, [_Mic(1)], coords)

            read.assert_not_called()


class _Closing(_Picking):
    def __init__(self, root, failed):
        super().__init__(root)
        self.micDict = OrderedDict([("mic_001", _Mic(1))])
        self.streamClosed = True
        self._failed = failed
        self.published = []

    def _getPublishedPickingMicIds(self):
        return set()

    def _getFinishedPickingMicNames(self):
        return {"mic_001"}

    def _updateOutputCoordSet(self, micList, streamMode):
        self.published.append([mic.getObjId() for mic in micList])
        if self._failed:
            self._recordFailedPicking(micList[0])

    def _updateStreamState(self, streamMode):
        pass

    def _getFirstJoinStep(self):
        return None

    def _streamingSleepOnWait(self):
        pass


class TestPickingAlwaysEnds(unittest.TestCase):
    def test_AMicrographWithoutCoordinatesIsReadOnce(self):
        with tempfile.TemporaryDirectory() as root:
            protocol = _Closing(root, failed=False)
            protocol.streamClosed = False

            ProtRelionAutopickBase._checkNewOutput(protocol)
            ProtRelionAutopickBase._checkNewOutput(protocol)

            self.assertEqual([[1]], protocol.published)

    def test_EveryMicrographFailedFailsTheRun(self):
        with tempfile.TemporaryDirectory() as root:
            with self.assertRaisesRegex(Exception, "Picking failed for every"):
                ProtRelionAutopickBase._checkNewOutput(_Closing(root, failed=True))

    def test_SomeMicrographsWithoutParticlesAreNoFailure(self):
        with tempfile.TemporaryDirectory() as root:
            protocol = _Closing(root, failed=False)

            ProtRelionAutopickBase._checkNewOutput(protocol)

            self.assertTrue(protocol.finished)


class TestFailedMicrographsArePickedAgain(unittest.TestCase):
    def test_AFailedMicrographLosesItsDoneMarker(self):
        with tempfile.TemporaryDirectory() as root:
            protocol = _Picking(root)
            marker = protocol._getExtraPath('DONE', 'mic_000001.TXT')
            os.makedirs(os.path.dirname(marker))
            open(marker, 'w').close()
            protocol._getMicDone = lambda mic: marker
            protocol.micDict["mic_001"] = _Mic(1)
            protocol._recordFailedPicking(_Mic(1))

            protocol._forgetDoneMarkersOfFailed(["mic_001"])

            self.assertFalse(os.path.exists(marker))

    def test_AMicrographPickedWellIsNoLongerFailed(self):
        with tempfile.TemporaryDirectory() as root:
            protocol = _Picking(root)
            protocol._recordFailedPicking(_Mic(1))

            protocol._forgetFailedPicking([1])

            self.assertEqual(set(), protocol._failedPickingMicIds())
            # Kept on disk: a later run reads the same list.
            self.assertEqual(set(), _Picking(root)._failedPickingMicIds())


if __name__ == '__main__':
    unittest.main()
