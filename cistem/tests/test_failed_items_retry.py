# **************************************************************************
# *
# * Movies and micrographs that failed in an earlier run are processed
# * again on Continue, and picking fails only when every micrograph failed.
# *
# **************************************************************************

import os
import tempfile
import unittest
from collections import OrderedDict

from cistem.protocols.protocol_picking import CistemProtFindParticles
from cistem.protocols.protocol_unblur import CistemProtUnblur

from .test_continue_after_interruption import _UnblurContinueHarness, _movieStep
from .test_streaming_architecture import _Mic, _StoredValue


class _UnblurRetryHarness(_UnblurContinueHarness):
    def __init__(self, resultExists, **kwargs):
        super().__init__(**kwargs)
        self._resultExists = resultExists

    def _movieResultExists(self, movie):
        return self._resultExists[movie.getObjId()]


class TestUnblurAlignsFailedMoviesAgain(unittest.TestCase):
    def test_AMovieWithoutAlignedResultIsProcessedAgain(self):
        # Both were processed by the earlier run; only movie 1 left its
        # aligned micrograph behind.
        protocol = _UnblurRetryHarness(
            {1: True, 2: False}, movieIds=[1, 2], publishedIds=[],
            prevSteps=[_movieStep(1, True), _movieStep(2, True)])

        CistemProtUnblur._restoreProcessedMoviesFromPersistentState(protocol)

        self.assertEqual({1, 2}, set(protocol._pendingMovies))
        self.assertIn(1, protocol.insertedDict)  # only published
        self.assertNotIn(2, protocol.insertedDict)  # gets a new step
        # Its old finished step must not count: it would be published
        # (and fail) before the new one runs.
        self.assertEqual({1}, CistemProtUnblur._getFinishedMovieIds(protocol))


class _PickStep:
    def __init__(self, micName, finished=True):
        self.funcName = _StoredValue("pickMicrographStep")
        self.argsStr = _StoredValue('["%s", {}]' % micName)
        self._finished = finished

    def isFinished(self):
        return self._finished


class _PickingHarness:
    def __init__(self, root, failedIds=()):
        self._root = root
        os.makedirs(os.path.join(root, 'extra'), exist_ok=True)
        self.micDict = OrderedDict()
        self._pendingMics = OrderedDict()
        self._prevSteps = [_PickStep("mic_001"), _PickStep("mic_002")]
        self._steps = []
        if failedIds:
            CistemProtFindParticles._writeFailedList(
                self, [_Mic(i, "mic_%03d" % i) for i in failedIds])

    def _getExtraPath(self, *paths):
        return os.path.join(self._root, 'extra', *paths)

    _getAllFailed = CistemProtFindParticles._getAllFailed
    _readFailedList = CistemProtFindParticles._readFailedList
    _writeFailedList = CistemProtFindParticles._writeFailedList
    _forgetFailedMics = CistemProtFindParticles._forgetFailedMics


class TestPickingPicksFailedMicrographsAgain(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self._tmp.cleanup()

    def test_AFailedMicrographIsPickedAgainAndForgottenAsFailed(self):
        protocol = _PickingHarness(self._tmp.name, failedIds=[2])
        ready = OrderedDict([("mic_001", _Mic(1, "mic_001")),
                             ("mic_002", _Mic(2, "mic_002"))])

        CistemProtFindParticles._takeAlreadyPickedMics(
            protocol, ready, {"mic_001", "mic_002"})

        self.assertEqual(["mic_001"], list(protocol.micDict))
        self.assertEqual(["mic_002"], list(ready))
        self.assertEqual({"mic_002"}, protocol._retriedMicNames)
        self.assertEqual([], CistemProtFindParticles._readFailedList(protocol))

    def test_EveryMicrographFailedFailsTheRun(self):
        protocol = _PickingHarness(self._tmp.name, failedIds=[1, 2])
        protocol._knownMicIds = {1, 2}

        with self.assertRaisesRegex(Exception, "Picking failed for every"):
            CistemProtFindParticles._failIfEveryMicrographFailed(protocol)

    def test_MicrographsWithoutParticlesAreNotFailures(self):
        protocol = _PickingHarness(self._tmp.name, failedIds=[1])
        protocol._knownMicIds = {1, 2}

        CistemProtFindParticles._failIfEveryMicrographFailed(protocol)


if __name__ == '__main__':
    unittest.main()
