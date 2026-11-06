"""Failure-contract regressions for Relion LoG auto-picking.

Unit tests intentionally do not run Relion or MPI. They validate the boundary
between a failed executable and a legitimately empty picking result.
"""

import os
import tempfile
from unittest import TestCase
from unittest.mock import patch

from relion.protocols.protocol_autopick_log import ProtRelionAutopickLoG


class _Mic:
    def __init__(self, obj_id):
        self._obj_id = obj_id

    def getObjId(self):
        return self._obj_id


class _LoGFailureHarness(ProtRelionAutopickLoG):
    """Exercise the real LoG batch method without a full Scipion project."""

    def __init__(self, root, failure=None, produced_ids=()):
        self.root = root
        self.failure = failure
        self.produced_ids = set(produced_ids)
        self.numberOfMpi = 2
        self.errors = []
        self.warnings = []
        self.jobs = []

    def _getExtraPath(self, *parts):
        os.makedirs(os.path.join(self.root, "extra"), exist_ok=True)
        return os.path.join(self.root, "extra", *parts)

    def _getTmpPath(self, *parts):
        return os.path.join(self.root, *parts)

    def _createTmpMicsDir(self, mic_list):
        path = self._getTmpPath("batch")
        os.makedirs(path, exist_ok=True)
        return path

    def runJob(self, program, args, **kwargs):
        self.jobs.append((program, args))
        if self.failure == "mpi":
            raise RuntimeError("2 total processes failed to start")

        # Even an empty STAR is a valid result: zero detected particles.
        for obj_id in self.produced_ids:
            path = os.path.join(kwargs["cwd"],
                                "mic_%06d_autopick.star" % obj_id)
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("data_\n\nloop_\n_rlnCoordinateX #1\n_rlnCoordinateY #2\n")

    def error(self, message):
        self.errors.append(message)

    def warning(self, message):
        self.warnings.append(message)


class TestRelionLoGFailureContract(TestCase):
    def _pick(self, protocol, mic_ids):
        # The STAR writer is unrelated to this contract and needs no input DB.
        with patch("relion.protocols.protocol_autopick.convert.createWriter"):
            protocol._pickMicrographList(
                [_Mic(obj_id) for obj_id in mic_ids],
                " --pickname autopick --LoG", 114, 138, 0, 2.5,
            )

    # A failed micrograph must never look like one without particles. Under
    # the streaming policy it is recorded as failed (not published, picked
    # again on Continue) and the run fails only when every one failed.

    def test_AnMpiLaunchFailureRecordsTheWholeBatchAsFailed(self):
        with tempfile.TemporaryDirectory() as root:
            protocol = _LoGFailureHarness(root, failure="mpi")

            self._pick(protocol, [19, 20, 21])

            self.assertEqual(1, len(protocol.jobs))
            self.assertEqual("relion_autopick_mpi", protocol.jobs[0][0])
            self.assertEqual({19, 20, 21}, protocol._failedPickingMicIds())
            self.assertTrue(protocol.errors)

    def test_ASuccessExitWithoutAnyStarRecordsTheMicrographsAsFailed(self):
        with tempfile.TemporaryDirectory() as root:
            protocol = _LoGFailureHarness(root)

            self._pick(protocol, [19, 20, 21])

            self.assertEqual({19, 20, 21}, protocol._failedPickingMicIds())

    def test_OnlyTheMicrographMissingFromABatchFails(self):
        with tempfile.TemporaryDirectory() as root:
            protocol = _LoGFailureHarness(root, produced_ids=(19, 21))

            self._pick(protocol, [19, 20, 21])

            self.assertEqual({20}, protocol._failedPickingMicIds())

    def test_ValidStarsWithZeroParticlesAreNotAnExecutionFailure(self):
        with tempfile.TemporaryDirectory() as root:
            protocol = _LoGFailureHarness(root, produced_ids=(19, 20, 21))

            self._pick(protocol, [19, 20, 21])

            for obj_id in (19, 20, 21):
                self.assertTrue(os.path.isfile(protocol._getTmpPath(
                    "mic_%06d_autopick.star" % obj_id
                )))
            self.assertEqual([], protocol.errors)
            self.assertEqual([], protocol.warnings)
            self.assertEqual(set(), protocol._failedPickingMicIds())
