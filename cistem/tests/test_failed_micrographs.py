# **************************************************************************
# *
# * A micrograph whose CTF could not be estimated is never published (no
# * -999 placeholder): Continue estimates it again, and the run fails only
# * when no micrograph at all could be estimated.
# *
# **************************************************************************

import unittest
from collections import OrderedDict

from pwem.objects import CTFModel

from cistem.protocols.protocol_ctffind import CistemProtCTFFind
from cistem.protocols.protocol_streaming_base import CistemStreamingBase

from .test_streaming_architecture import (
    _CtffindOutputCheckHarness,
    _CtffindResumeHarness,
    _LogicalMicrographSet,
    _Mic,
    _StoredValue,
)


class _Step:
    def __init__(self, micName, finished=True):
        self.funcName = _StoredValue("estimateCtfStep")
        self.argsStr = _StoredValue('["%s"]' % micName)
        self._finished = finished

    def isFinished(self):
        return self._finished


def _ctf(defocusU, defocusV):
    ctf = CTFModel()
    ctf.setStandardDefocus(defocusU, defocusV, 45.0)
    return ctf


class _Program:
    def __init__(self, ctf):
        self.ctf = ctf

    def parseOutputAsCtf(self, ctfFn, rotAvgFn=None, psdFile=None):
        return self.ctf


class _ModelHarness:
    def __init__(self, ctf):
        self._ctfProgram = _Program(ctf)

    def _getPsdPath(self, mic):
        return "psd.mrc"

    def _getCtfOutPath(self, mic):
        return "ctf.txt"

    def _getCtfAvrotPath(self, mic):
        return "ctf_avrot.txt"


class TestAFailedEstimationIsNeverPublished(unittest.TestCase):
    def test_AnUnusableResultIsRefusedInsteadOfPublishedAsMinus999(self):
        failed = CTFModel()
        failed.setWrongDefocus()
        protocol = _ModelHarness(failed)

        with self.assertRaises(ValueError):
            CistemProtCTFFind._createCtfModel(protocol, _Mic(1, "mic_001"))

    def test_AGoodResultIsPublishedAndItsFileKept(self):
        protocol = _ModelHarness(_ctf(12000.0, 11000.0))

        ctf = CistemProtCTFFind._createCtfModel(protocol, _Mic(1, "mic_001"))

        self.assertEqual(12000.0, ctf.getDefocusU())


class _UsabilityResumeHarness(_CtffindResumeHarness):
    def __init__(self, inputSet, usable):
        super().__init__(inputSet)
        self._usable = usable

    def _ctfResultIsUsable(self, mic):
        return self._usable[mic.getMicName()]


class TestContinueEstimatesFailedMicrographsAgain(unittest.TestCase):
    def test_OnlyMicrographsWithAnUnusableResultAreEstimatedAgain(self):
        mic1, mic2 = _Mic(1, "mic_001"), _Mic(2, "mic_002")
        protocol = _UsabilityResumeHarness(
            _LogicalMicrographSet([mic1, mic2], streamClosed=True),
            usable={"mic_001": False, "mic_002": True})
        protocol._prevSteps = [_Step("mic_001"), _Step("mic_002")]

        CistemProtCTFFind.stepsGeneratorStep(protocol)

        self.assertEqual(["mic_001"], protocol.scheduledMicNames)
        self.assertIn("mic_002", protocol.micDict)

    def test_ARetriedMicrographWaitsForItsNewStep(self):
        # Its step from the earlier run finished with an unusable result:
        # publishing on that one would drop the micrograph again.
        protocol = _CtffindOutputCheckHarness()
        protocol._prevSteps, protocol._steps = protocol._steps, []
        protocol._retriedMicNames = {"mic_001"}

        CistemProtCTFFind._checkNewOutput(protocol)
        self.assertEqual([], protocol.publishCalls)

        protocol._steps = [_Step("mic_001")]
        CistemProtCTFFind._checkNewOutput(protocol)
        self.assertEqual([["mic_001"]], protocol.publishCalls)


class _Output:
    def __init__(self, size):
        self._size = size

    def getSize(self):
        return self._size


class _ClosingHarness(_CtffindOutputCheckHarness):
    def __init__(self, publishedCount, knownIds):
        super().__init__()
        self.streamClosed = True
        self._knownMicIds = set(knownIds)
        self.outputCTF = _Output(publishedCount)

    def _updateStreamState(self, streamMode):
        pass

    def _failIfNoCtfWasEstimated(self):
        return CistemProtCTFFind._failIfNoCtfWasEstimated(self)


class TestTheRunFailsOnlyWhenEveryMicrographFailed(unittest.TestCase):
    def test_NothingPublishedAtTheEndFailsTheRun(self):
        with self.assertRaisesRegex(Exception, "No CTF could be estimated"):
            CistemProtCTFFind._checkNewOutput(_ClosingHarness(0, {1}))

    def test_SomeFailuresDoNotFailTheRun(self):
        CistemProtCTFFind._checkNewOutput(_ClosingHarness(1, {1, 2}))

    def test_AnEmptyStreamIsNotAFailure(self):
        protocol = _ClosingHarness(0, set())
        protocol.micDict = OrderedDict()
        protocol._steps = []

        CistemProtCTFFind._checkNewOutput(protocol)

        self.assertTrue(protocol.finished)


class TestStepsOfThisRunOnly(unittest.TestCase):
    def test_PreviousRunStepsCanBeLeftOut(self):
        protocol = CistemStreamingBase()
        protocol._prevSteps = [_Step("mic_001")]
        protocol._steps = [_Step("mic_002")]

        self.assertEqual({"mic_002"}, protocol._collectStepArgKeys(
            ("estimateCtfStep",), includePrevious=False))
        self.assertEqual({"mic_001", "mic_002"},
                         protocol._collectStepArgKeys(("estimateCtfStep",)))


if __name__ == '__main__':
    unittest.main()
