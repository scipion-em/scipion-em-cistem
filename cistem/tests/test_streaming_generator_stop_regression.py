# **************************************************************************
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
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""The streaming generator must abandon its loop once the run is over.

A failed step makes pyworkflow mark the protocol as FAILED and the step
executor break out of its own loop - and then join every thread it
started, the steps generator among them. A generator that keeps polling
is never joined, so the run hangs for good with nothing left to do.

All three cisTEM streaming protocols run a polling loop of their own, so
each one is checked here.
"""
import unittest

import pyworkflow.protocol.constants as cons

from cistem.protocols.protocol_ctffind import CistemProtCTFFind
from cistem.protocols.protocol_picking import CistemProtFindParticles
from cistem.protocols.protocol_streaming_base import CistemStreamingBase
from cistem.protocols.protocol_unblur import CistemProtUnblur


# A real hang cannot be asserted on, so every harness caps its own
# polling: reaching the cap is the failure signal.
MAX_POLLS = 50


class _Status:
    """Mimics the String attribute pyworkflow keeps the status in."""

    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _GeneratorHarness:
    """Collaborators every cisTEM generator leans on, stubbed out."""

    _streamingMustStop = CistemStreamingBase._streamingMustStop

    def __init__(self, status=cons.STATUS_RUNNING, closeAfter=None):
        self.status = _Status(status)
        self.finished = False
        self.polls = 0
        self.sleeps = 0
        self.insertedFinalSteps = 0
        self.insertedOutputSteps = 0
        self._closeAfter = closeAfter

    # -- the generator's own bookkeeping ---------------------------------
    def _defineCtfParamsDict(self):
        pass

    def _insertInitialSteps(self):
        return []

    def _restoreProcessedMicsFromPersistentState(self):
        pass

    _restoreProcessedMoviesFromPersistentState = (
        _restoreProcessedMicsFromPersistentState
    )

    # -- the polling body ------------------------------------------------
    def _checkNewInput(self):
        self.polls += 1
        if self.polls > MAX_POLLS:
            raise AssertionError(
                "The generator polled %d times without leaving its loop: "
                "the executor's join() would never return." % self.polls
            )

    def _checkNewOutput(self):
        if self._closeAfter is not None and self.polls >= self._closeAfter:
            self.finished = True

    def _getStreamingSleepOnWait(self):
        return 1

    def _streamingSleepOnWait(self):
        self.sleeps += 1

    # -- unblur's terminal scheduling ------------------------------------
    def _insertFinalSteps(self, deps):
        self.insertedFinalSteps += 1
        return []

    def _insertFunctionStep(self, *args, **kwargs):
        self.insertedOutputSteps += 1
        return 1


class _StopsWhenTheRunIsOver:
    """Shared expectations, run against each protocol's own generator."""

    protocolClass = None

    def _runGenerator(self, harness):
        self.protocolClass.stepsGeneratorStep(harness)

    def testFailedRunLeavesTheLoopImmediately(self):
        harness = _GeneratorHarness(status=cons.STATUS_FAILED)

        self._runGenerator(harness)

        self.assertEqual(
            harness.polls,
            0,
            "A FAILED protocol must not poll even once: the executor is "
            "already joining this very thread.",
        )

    def testAbortedRunLeavesTheLoopImmediately(self):
        harness = _GeneratorHarness(status=cons.STATUS_ABORTED)

        self._runGenerator(harness)

        self.assertEqual(harness.polls, 0)

    def testFailureMidStreamStopsTheLoop(self):
        harness = _GeneratorHarness()
        realCheck = harness._checkNewInput

        def failOnThirdPoll():
            realCheck()
            if harness.polls == 3:
                harness.status = _Status(cons.STATUS_FAILED)

        harness._checkNewInput = failOnThirdPoll

        self._runGenerator(harness)

        self.assertEqual(
            harness.polls,
            3,
            "The loop must stop on the poll right after the failure, not "
            "keep spinning until the producer closes the stream.",
        )

    def testHealthyRunStillPollsUntilTheStreamCloses(self):
        harness = _GeneratorHarness(closeAfter=4)

        self._runGenerator(harness)

        self.assertEqual(
            harness.polls,
            4,
            "A running protocol must keep polling: the stop guard must "
            "not short-circuit a healthy stream.",
        )

    def testAPlainStringStatusIsUnderstood(self):
        """status is not always wrapped: accept the bare value too."""
        harness = _GeneratorHarness()
        harness.status = cons.STATUS_FAILED

        self._runGenerator(harness)

        self.assertEqual(harness.polls, 0)

    def testAMissingStatusDoesNotBreakTheLoop(self):
        harness = _GeneratorHarness(closeAfter=1)
        del harness.status

        self._runGenerator(harness)

        self.assertEqual(harness.polls, 1)


class TestCtffindGeneratorStops(_StopsWhenTheRunIsOver, unittest.TestCase):
    protocolClass = CistemProtCTFFind


class TestPickingGeneratorStops(_StopsWhenTheRunIsOver, unittest.TestCase):
    protocolClass = CistemProtFindParticles


class TestUnblurGeneratorStops(_StopsWhenTheRunIsOver, unittest.TestCase):
    protocolClass = CistemProtUnblur

    def testAFailedRunSchedulesNoTerminalSteps(self):
        """Unblur inserts its output step after the loop.

        A step inserted once the run has failed is one the executor will
        never get to, so the generator must not reach that code at all.
        """
        harness = _GeneratorHarness(status=cons.STATUS_FAILED)

        self._runGenerator(harness)

        self.assertEqual(harness.insertedFinalSteps, 0)
        self.assertEqual(
            harness.insertedOutputSteps,
            0,
            "Scheduling createOutputStep on a failed run only adds a step "
            "that can never run.",
        )

    def testAHealthyRunStillSchedulesItsTerminalSteps(self):
        harness = _GeneratorHarness(closeAfter=2)

        self._runGenerator(harness)

        self.assertEqual(harness.insertedFinalSteps, 1)
        self.assertEqual(harness.insertedOutputSteps, 1)


if __name__ == '__main__':
    unittest.main()


class _StuckClosedSet:
    """A producer that closed declaring more items than it ever shows.

    The declared size never comes down and the missing row never turns
    up: the view is terminally inconsistent, for good.
    """

    def __init__(self, declaredSize=100, visibleIds=None):
        self._declaredSize = declaredSize
        self._visibleIds = list(visibleIds if visibleIds is not None
                                else range(1, 100))

    def getSize(self):
        return self._declaredSize

    def isStreamClosed(self):
        return True

    def getUniqueValues(self, attributes, where=None):
        return list(self._visibleIds) if where is None else []

    def close(self):
        pass


class _TerminalStallHarness(CistemStreamingBase):
    """Reconciles a closed-but-inconsistent stream, poll after poll."""

    def __init__(self, activeWork=False):
        self._lastInputId = 0
        self._activeWork = activeWork

    def _hasActiveStreamingWork(self):
        return self._activeWork

    def poll(self, inputSet):
        return self._reconcileClosedStreamIds(inputSet, [], set(), True, '_lastInputId')


class TestTerminalInconsistencyDoesNotHangForever(unittest.TestCase):
    """A closed producer whose view never becomes consistent must not
    leave the protocol polling for the rest of time."""

    def _pollUntilRaises(self, harness, inputSet, limit=200):
        for poll in range(limit):
            try:
                harness.poll(inputSet)
            except RuntimeError as error:
                return poll + 1, str(error)

        return None, None

    def testAPermanentlyInconsistentViewEventuallyFails(self):
        polls, message = self._pollUntilRaises(_TerminalStallHarness(),
                                               _StuckClosedSet())

        self.assertIsNotNone(
            polls,
            "The producer closed declaring 100 items and only 99 are ever "
            "visible: polling for that hundredth row never ends.",
        )
        self.assertIn('99', message)
        self.assertIn('100', message)

    def testItGivesTheViewSeveralChancesFirst(self):
        """A lagging view usually catches up; do not fail on poll one."""
        polls, _ = self._pollUntilRaises(_TerminalStallHarness(),
                                         _StuckClosedSet())

        self.assertGreater(
            polls,
            3,
            "Giving up almost immediately would turn an ordinary lag into "
            "a failed protocol.",
        )

    def testProgressResetsTheCount(self):
        harness = _TerminalStallHarness()
        inputSet = _StuckClosedSet()

        for _ in range(5):
            harness.poll(inputSet)

        inputSet._visibleIds.append(100)
        harness.poll(inputSet)

        self.assertEqual(
            harness._terminalStallCount,
            0,
            "The view did become consistent; nothing is stalled.",
        )

    def testWorkInFlightIsAlsoProgress(self):
        """A long scientific round must not be mistaken for a stall."""
        polls, _ = self._pollUntilRaises(
            _TerminalStallHarness(activeWork=True), _StuckClosedSet(),
            limit=60)

        self.assertIsNone(
            polls,
            "Work was in flight the whole time: that is progress, however "
            "long it takes.",
        )

    def testAConsistentViewIsNeverAffected(self):
        harness = _TerminalStallHarness()
        inputSet = _StuckClosedSet(declaredSize=99)

        for _ in range(50):
            _, consistent = harness.poll(inputSet)

        self.assertTrue(consistent)
