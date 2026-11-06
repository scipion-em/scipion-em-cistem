# **************************************************************************
# *
# * Continue after an interrupted run: every item ends published or
# * processed again, and a stream that never completes is reported.
# *
# **************************************************************************

import unittest

from cistem.protocols.protocol_ctffind import CistemProtCTFFind
from cistem.protocols.protocol_picking import CistemProtFindParticles
from cistem.protocols.protocol_streaming_base import CistemStreamingBase
from cistem.protocols.protocol_unblur import CistemProtUnblur

from .logical_set_fakes import LogicalSetFake
from .test_streaming_architecture import (
    _CtffindResumeHarness,
    _LogicalMicrographSet,
    _LogicalMovie,
    _LogicalMovieSet,
    _Mic,
    _StoredValue,
    _UnblurPointer,
)


class _Step:
    """A persisted step as pyworkflow restores it on Continue."""

    def __init__(self, funcName, argsStr, finished):
        self.funcName = _StoredValue(funcName)
        self.argsStr = _StoredValue(argsStr)
        self._finished = finished

    def isFinished(self):
        return self._finished


def _ctfStep(micName, finished):
    return _Step("estimateCtfStep", '["%s"]' % micName, finished)


def _movieStep(movieId, finished):
    return _Step("processMovieStep", '[{"object.id": %d}, false]' % movieId,
                 finished)


class TestStepsOfAnEarlierRunThatDidNotFinish(unittest.TestCase):
    """pyworkflow only reruns the steps _insertAllSteps creates: the ones a
    generator inserted in an earlier run stay in _prevSteps and never run
    again, so an unfinished one there is not "scheduled" any more."""

    def test_OnlyFinishedStepsOfAnEarlierRunCountAsScheduled(self):
        protocol = CistemStreamingBase()
        protocol._prevSteps = [_ctfStep("mic_001", finished=False),
                               _ctfStep("mic_002", finished=True)]
        protocol._steps = [_ctfStep("mic_003", finished=False)]

        scheduled = protocol._collectStepArgKeys(("estimateCtfStep",),
                                                 onlyFinished=False)

        # mic_003 belongs to this run and will still run.
        self.assertEqual({"mic_002", "mic_003"}, scheduled)

    def test_CtffindProcessesAgainAMicrographWhoseStepWasInterrupted(self):
        mic1, mic2 = _Mic(1, "mic_001"), _Mic(2, "mic_002")
        protocol = _CtffindResumeHarness(
            _LogicalMicrographSet([mic1, mic2], streamClosed=True))
        protocol._prevSteps = [_ctfStep("mic_001", finished=False),
                               _ctfStep("mic_002", finished=True)]

        CistemProtCTFFind.stepsGeneratorStep(protocol)

        # mic_001 is estimated again; mic_002 only waits to be published.
        self.assertEqual(["mic_001"], protocol.scheduledMicNames)
        self.assertIn("mic_002", protocol.micDict)


class _UnblurContinueHarness(CistemStreamingBase):
    def __init__(self, movieIds, publishedIds, prevSteps):
        self.inputMovies = _UnblurPointer(_LogicalMovieSet(
            [_LogicalMovie(i) for i in movieIds], streamClosed=True))
        self.outputMicrographs = _LogicalMovieSet(
            [_LogicalMovie(i) for i in publishedIds], streamClosed=False)
        self._prevSteps = prevSteps
        self._steps = []
        self.insertedDict = {}
        self.listOfMovies = []

    def _getPublishedMovieIds(self):
        return CistemProtUnblur._getPublishedMovieIds(self)

    def _getScheduledMovieIds(self):
        return CistemProtUnblur._getScheduledMovieIds(self)

    def _movieResultExists(self, movie):
        # What an earlier run aligned is there unless a test says not.
        return True

    def debug(self, *args, **kwargs):
        pass


class TestUnblurLosesNoMovieOnContinue(unittest.TestCase):
    def test_EveryMovieEndsPublishedOrProcessedAgain(self):
        # 1: aligned but never published; 2: its step was interrupted;
        # 3: already published.
        protocol = _UnblurContinueHarness(
            movieIds=[1, 2, 3], publishedIds=[3],
            prevSteps=[_movieStep(1, finished=True),
                       _movieStep(2, finished=False)])

        CistemProtUnblur._restoreProcessedMoviesFromPersistentState(protocol)
        CistemProtUnblur._loadInputList(protocol)

        self.assertEqual({1, 2}, set(protocol._pendingMovies))
        # 1 only has to be published; 2 gets a new step.
        self.assertIn(1, protocol.insertedDict)
        self.assertNotIn(2, protocol.insertedDict)


class _StalledStream(LogicalSetFake):
    """A closed producer that declares more items than it ever shows."""

    def __init__(self, visibleIds, declared):
        super().__init__([_Mic(i, "mic_%03d" % i) for i in visibleIds],
                         streamClosed=True)
        self._declared = declared

    def getSize(self):
        return self._declared


class _IdleProtocol(CistemStreamingBase):
    def __init__(self, protocolClass, **work):
        self._protocolClass = protocolClass
        for name, value in work.items():
            setattr(self, name, value)

    def _hasActiveStreamingWork(self):
        return self._protocolClass._hasActiveStreamingWork(self)


class TestAStreamThatNeverCompletesIsReported(unittest.TestCase):
    def _poll(self, protocol, stream, attr, times):
        for _ in range(times):
            protocol._reconcileClosedStreamIds(stream, [], {1}, True, attr)

    def test_CtffindWithNothingInFlightFailsInsteadOfWaitingForever(self):
        protocol = _IdleProtocol(CistemProtCTFFind, micDict={},
                                 _pendingMics={})
        stream = _StalledStream([1], declared=3)

        with self.assertRaisesRegex(RuntimeError, "declaring 3 items"):
            self._poll(protocol, stream, '_lastInputId',
                       CistemStreamingBase.TERMINAL_STALL_POLLS)

    def test_WorkInFlightIsNeverCountedAsAStall(self):
        protocol = _IdleProtocol(CistemProtCTFFind,
                                 micDict={"mic_001": _Mic(1, "mic_001")},
                                 _pendingMics={})
        stream = _StalledStream([1], declared=3)

        self._poll(protocol, stream, '_lastInputId',
                   CistemStreamingBase.TERMINAL_STALL_POLLS * 3)

    def test_EveryStreamingProtocolKnowsWhenItIsIdle(self):
        for protocolClass, work in (
                (CistemProtCTFFind, dict(micDict={}, _pendingMics={})),
                (CistemProtFindParticles, dict(micDict={}, _pendingMics={})),
                (CistemProtUnblur, dict(_pendingMovies={}))):
            idle = _IdleProtocol(protocolClass, **work)
            self.assertFalse(idle._hasActiveStreamingWork(), protocolClass)

    def test_PickingCountsItsTwoStreamsApart(self):
        # Micrographs and CTFs are polled alternately; one shared counter
        # was reset by the other stream on every poll and never fired.
        protocol = _IdleProtocol(CistemProtFindParticles, micDict={},
                                 _pendingMics={})
        mics = _StalledStream([1], declared=3)
        ctfs = _StalledStream([1], declared=3)

        with self.assertRaisesRegex(RuntimeError, "declaring 3 items"):
            for _ in range(CistemStreamingBase.TERMINAL_STALL_POLLS):
                self._poll(protocol, mics, '_lastMicId', 1)
                self._poll(protocol, ctfs, '_lastCtfId', 1)


if __name__ == '__main__':
    unittest.main()
