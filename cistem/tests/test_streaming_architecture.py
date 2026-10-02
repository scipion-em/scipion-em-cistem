# **************************************************************************
# *
# * Focused regression tests for cisTEM streaming architecture.
# *
# **************************************************************************

import unittest

from cistem.protocols.protocol_ctffind import (
    CistemProtCTFFind,
)


class _Mic:
    def __init__(self, objId, micName):
        self._objId = objId
        self._micName = micName

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return self._micName

    def clone(self):
        return _Mic(
            self._objId,
            self._micName,
        )


class _LogicalMicrographSet:
    def __init__(self, items, streamClosed=False):
        self._items = list(items)
        self._streamClosed = streamClosed
        self.loadCalls = 0

    def getFileName(self):
        raise AssertionError(
            "Streaming discovery must not depend on "
            "a SQLite/storage filename."
        )

    def loadAllProperties(self):
        self.loadCalls += 1

    def iterItems(self):
        return iter(self._items)

    def isStreamClosed(self):
        return self._streamClosed


class _CtffindStreamingHarness:
    def __init__(self, inputSet):
        self._inputSet = inputSet
        self.micDict = {}
        self.debugMessages = []

    def getInputMicrographs(self):
        return self._inputSet

    def debug(self, message):
        self.debugMessages.append(message)

    def _loadSet(
            self,
            inputSet,
            setClass,
            getKeyFunc,
    ):
        return CistemProtCTFFind._loadSet(
            self,
            inputSet,
            setClass,
            getKeyFunc,
        )


class TestCistemCtffindStreamingArchitecture(
        unittest.TestCase,
):
    def test_CtffindStreamingLoadsLogicalMicrographsWithoutStorageFilename(
            self,
    ):
        inputSet = _LogicalMicrographSet(
            [
                _Mic(1, "mic_001"),
                _Mic(2, "mic_002"),
            ],
            streamClosed=False,
        )

        protocol = _CtffindStreamingHarness(
            inputSet
        )

        newMics, streamClosed = (
            CistemProtCTFFind
            ._loadInputList(
                protocol
            )
        )

        self.assertEqual(
            list(newMics),
            [
                "mic_001",
                "mic_002",
            ],
        )

        self.assertFalse(
            streamClosed
        )

        self.assertEqual(
            inputSet.loadCalls,
            1,
        )


class _JoinStep:
    def __init__(self):
        self.prerequisites = []

    def addPrerequisites(self, *deps):
        self.prerequisites.extend(deps)


class _CtffindInputCheckHarness(
        _CtffindStreamingHarness
):
    def __init__(self, inputSet):
        super().__init__(
            inputSet
        )
        self.streamClosed = False
        self.joinStep = _JoinStep()
        self.insertedMicNames = []
        self.updateStepsCalls = 0

    def _loadInputList(self):
        return (
            CistemProtCTFFind
            ._loadInputList(
                self
            )
        )

    def _insertNewMicsSteps(
            self,
            newMics,
    ):
        newMics = list(newMics)

        self.insertedMicNames.extend(
            mic.getMicName()
            for mic in newMics
        )

        return [
            101 + index
            for index, _ in enumerate(
                newMics
            )
        ]

    def _getFirstJoinStep(self):
        return self.joinStep

    def updateSteps(self):
        self.updateStepsCalls += 1


class TestCistemCtffindStreamingInputChecks(
        unittest.TestCase,
):
    def test_CtffindStreamingChecksLogicalInputWithoutFilesystemMtime(
            self,
    ):
        inputSet = _LogicalMicrographSet(
            [
                _Mic(1, "mic_001"),
                _Mic(2, "mic_002"),
            ],
            streamClosed=False,
        )

        protocol = (
            _CtffindInputCheckHarness(
                inputSet
            )
        )

        CistemProtCTFFind._checkNewInput(
            protocol
        )

        self.assertEqual(
            protocol.insertedMicNames,
            [
                "mic_001",
                "mic_002",
            ],
        )

        self.assertEqual(
            protocol.joinStep.prerequisites,
            [
                101,
                102,
            ],
        )

        self.assertEqual(
            protocol.updateStepsCalls,
            1,
        )

        self.assertFalse(
            protocol.streamClosed
        )


class _StoredValue:
    def __init__(self, value):
        self._value = value

    def get(self, default=None):
        if self._value is None:
            return default
        return self._value


class _FinishedCtfStep:
    def __init__(self, micName):
        self.funcName = _StoredValue(
            "estimateCtfStep"
        )
        self.argsStr = _StoredValue(
            '["%s"]' % micName
        )

    def isFinished(self):
        return True


class _CtffindOutputCheckHarness:
    def __init__(self):
        mic = _Mic(
            1,
            "mic_001",
        )

        self.micDict = {
            mic.getMicName(): mic,
        }
        self.streamClosed = False
        self.finished = False
        self._steps = [
            _FinishedCtfStep(
                mic.getMicName()
            ),
        ]

        self.publishedMicNames = set()
        self.publishCalls = []
        self.sleepCalls = 0

    def _readDoneList(self):
        raise AssertionError(
            "CTFFind streaming output must not depend on "
            "DONE/all.TXT."
        )

    def _isMicDone(self, mic):
        raise AssertionError(
            "CTFFind streaming output must not depend on "
            "DONE/mic_*.TXT."
        )

    def _getPublishedCtfMicNames(self):
        return set(
            self.publishedMicNames
        )

    def _getFinishedCtfMicNames(self):
        return (
            CistemProtCTFFind
            ._getFinishedCtfMicNames(
                self
            )
        )

    def _updateOutputCTFSet(
            self,
            micList,
            streamMode,
    ):
        micList = list(
            micList
        )

        names = [
            mic.getMicName()
            for mic in micList
        ]

        self.publishCalls.append(
            names
        )

        self.publishedMicNames.update(
            names
        )

        return micList

    def _getFirstJoinStep(self):
        return None

    def _streamingSleepOnWait(self):
        self.sleepCalls += 1

    def debug(self, message):
        pass


class TestCistemCtffindStreamingCompletion(
        unittest.TestCase,
):
    def test_CtffindPublishesFinishedStepWithoutDoneSidecars(
            self,
    ):
        protocol = (
            _CtffindOutputCheckHarness()
        )

        CistemProtCTFFind._checkNewOutput(
            protocol
        )

        CistemProtCTFFind._checkNewOutput(
            protocol
        )

        self.assertEqual(
            protocol.publishCalls,
            [
                [
                    "mic_001",
                ],
            ],
        )

        self.assertEqual(
            protocol.sleepCalls,
            1,
        )


class _CtffindProcessingHarness:
    def __init__(self):
        mic = _Mic(
            1,
            "mic_001",
        )

        self.micDict = {
            mic.getMicName(): mic,
        }
        self.estimatedMicNames = []
        self.batchEstimatedMicNames = []
        self.initialIds = []

    def isContinued(self):
        return True

    def _getMicrographDone(self, mic):
        raise AssertionError(
            "CTFFind processing must not use DONE/mic_*.TXT."
        )

    def _estimateCTF(self, mic, *args):
        self.estimatedMicNames.append(
            mic.getMicName()
        )

    def _estimateCtfList(self, micList, *args):
        self.batchEstimatedMicNames.extend(
            mic.getMicName()
            for mic in micList
        )

    def info(self, message):
        pass


class TestCistemCtffindStreamingProcessing(
        unittest.TestCase,
):
    def test_CtffindSingleProcessingDoesNotWriteDoneSidecar(
            self,
    ):
        protocol = (
            _CtffindProcessingHarness()
        )

        CistemProtCTFFind.estimateCtfStep(
            protocol,
            "mic_001",
        )

        self.assertEqual(
            protocol.estimatedMicNames,
            [
                "mic_001",
            ],
        )

    def test_CtffindBatchProcessingDoesNotWriteDoneSidecars(
            self,
    ):
        protocol = (
            _CtffindProcessingHarness()
        )

        CistemProtCTFFind.estimateCtfListStep(
            protocol,
            [
                "mic_001",
            ],
        )

        self.assertEqual(
            protocol.batchEstimatedMicNames,
            [
                "mic_001",
            ],
        )


class _CtffindInitialStepsHarness:
    def _getExtraPath(self, *parts):
        raise AssertionError(
            "CTFFind streaming initialization must not create "
            "or resolve a DONE sidecar directory."
        )


class TestCistemCtffindStreamingInitialization(
        unittest.TestCase,
):
    def test_CtffindInitialStepsDoNotCreateDoneDirectory(
            self,
    ):
        protocol = (
            _CtffindInitialStepsHarness()
        )

        result = (
            CistemProtCTFFind
            ._insertInitialSteps(
                protocol
            )
        )

        self.assertEqual(
            result,
            [],
        )


class _FalseValue:
    def __bool__(self):
        return False


class _CtffindInsertHarness:
    recalculate = _FalseValue()

    def __init__(self):
        self.insertCalls = []
        self.legacyInitialCalls = 0
        self.legacyLoadCalls = 0
        self.legacyFinalCalls = 0

    def _defineCtfParamsDict(self):
        pass

    def resumableStepGeneratorStep(
            self,
            timestamp,
    ):
        pass

    def _insertInitialSteps(self):
        self.legacyInitialCalls += 1
        return []

    def _loadInputList(self):
        self.legacyLoadCalls += 1
        return {}, False

    def _insertNewMicsSteps(self, mics):
        return []

    def _insertFinalSteps(self, deps):
        self.legacyFinalCalls += 1
        return deps

    def _getFirstJoinStepName(self):
        return "createOutputStep"

    def _insertFunctionStep(
            self,
            funcName,
            *args,
            **kwargs,
    ):
        if callable(funcName):
            funcName = funcName.__name__

        self.insertCalls.append(
            {
                "funcName": funcName,
                "args": args,
                "kwargs": kwargs,
            }
        )

        return len(
            self.insertCalls
        )


class TestCistemCtffindGeneratorOrchestration(
        unittest.TestCase,
):
    def test_CtffindInsertAllStepsUsesSingleResumableGenerator(
            self,
    ):
        protocol = _CtffindInsertHarness()

        CistemProtCTFFind._insertAllSteps(
            protocol
        )

        self.assertEqual(
            protocol.legacyInitialCalls,
            0,
        )
        self.assertEqual(
            protocol.legacyLoadCalls,
            0,
        )
        self.assertEqual(
            protocol.legacyFinalCalls,
            0,
        )

        self.assertEqual(
            len(protocol.insertCalls),
            1,
        )

        self.assertEqual(
            protocol.insertCalls[0]["funcName"],
            "resumableStepGeneratorStep",
        )

        self.assertFalse(
            protocol.insertCalls[0]["kwargs"].get(
                "needsGPU",
                True,
            )
        )


class _PublishedCtf:
    def __init__(self, mic):
        self._mic = mic

    def getMicrograph(self):
        return self._mic


class _LogicalCtfOutput:
    def __init__(self, ctfs):
        self._ctfs = list(ctfs)

    def iterItems(self):
        return iter(self._ctfs)


class _CtffindResumeHarness:
    def _insertInitialSteps(self):
        return []

    def __init__(
            self,
            inputSet,
            publishedMics=None,
            existingSteps=None,
    ):
        self._inputSet = inputSet
        self.outputCTF = _LogicalCtfOutput(
            [
                _PublishedCtf(mic)
                for mic in (
                    publishedMics
                    or []
                )
            ]
        )
        self._steps = list(
            existingSteps
            or []
        )
        self.scheduledMicNames = []

    def _defineCtfParamsDict(self):
        pass

    def _getPublishedCtfMicNames(self):
        return (
            CistemProtCTFFind
            ._getPublishedCtfMicNames(
                self
            )
        )

    def _getFinishedCtfMicNames(self):
        return (
            CistemProtCTFFind
            ._getFinishedCtfMicNames(
                self
            )
        )

    def _restoreProcessedMicsFromPersistentState(
            self,
    ):
        return (
            CistemProtCTFFind
            ._restoreProcessedMicsFromPersistentState(
                self
            )
        )

    def getInputMicrographs(self):
        return self._inputSet

    def _loadSet(
            self,
            inputSet,
            setClass,
            getKeyFunc,
    ):
        return CistemProtCTFFind._loadSet(
            self,
            inputSet,
            setClass,
            getKeyFunc,
        )

    def _loadInputList(self):
        return (
            CistemProtCTFFind
            ._loadInputList(
                self
            )
        )

    def _getFirstJoinStep(self):
        return None

    def _insertNewMicsSteps(
            self,
            newMics,
    ):
        newMics = list(newMics)

        self.scheduledMicNames.extend(
            mic.getMicName()
            for mic in newMics
        )

        for mic in newMics:
            self.micDict[
                mic.getMicName()
            ] = mic

        return []

    def updateSteps(self):
        pass

    def _checkNewInput(self):
        CistemProtCTFFind._checkNewInput(
            self
        )

    def _checkNewOutput(self):
        self.finished = True

    def _getStreamingSleepOnWait(self):
        return 0

    def _streamingSleepOnWait(self):
        pass

    def debug(self, message):
        pass


class TestCistemCtffindGeneratorResumeSafety(
        unittest.TestCase,
):
    def test_CtffindGeneratorResumeDoesNotReschedulePublishedMicrographs(
            self,
    ):
        mic1 = _Mic(
            1,
            "mic_001",
        )
        mic2 = _Mic(
            2,
            "mic_002",
        )

        protocol = _CtffindResumeHarness(
            _LogicalMicrographSet(
                [
                    mic1,
                    mic2,
                ],
                streamClosed=True,
            ),
            publishedMics=[
                mic1,
            ],
        )

        CistemProtCTFFind.stepsGeneratorStep(
            protocol
        )

        self.assertEqual(
            protocol.scheduledMicNames,
            [
                "mic_002",
            ],
        )

    def test_CtffindGeneratorResumeDoesNotRescheduleFinishedUnpublishedSteps(
            self,
    ):
        mic1 = _Mic(
            1,
            "mic_001",
        )
        mic2 = _Mic(
            2,
            "mic_002",
        )

        protocol = _CtffindResumeHarness(
            _LogicalMicrographSet(
                [
                    mic1,
                    mic2,
                ],
                streamClosed=True,
            ),
            existingSteps=[
                _FinishedCtfStep(
                    "mic_001"
                ),
            ],
        )

        CistemProtCTFFind.stepsGeneratorStep(
            protocol
        )

        self.assertEqual(
            protocol.scheduledMicNames,
            [
                "mic_002",
            ],
        )


class _CtffindResumePublishHarness(
        _CtffindResumeHarness
):
    def __init__(
            self,
            inputSet,
            existingSteps=None,
    ):
        super().__init__(
            inputSet,
            publishedMics=[],
            existingSteps=existingSteps,
        )
        self.publishCalls = []

    def _updateOutputCTFSet(
            self,
            micList,
            streamMode,
    ):
        micList = list(
            micList
        )

        self.publishCalls.append(
            [
                mic.getMicName()
                for mic in micList
            ]
        )

        return micList

    def _updateStreamState(
            self,
            streamMode,
    ):
        pass

    def _getPublishedCtfMicNames(self):
        published = set()

        for batch in self.publishCalls:
            published.update(
                batch
            )

        return published

    def _checkNewOutput(self):
        return (
            CistemProtCTFFind
            ._checkNewOutput(
                self
            )
        )


class TestCistemCtffindGeneratorResumePublication(
        unittest.TestCase,
):
    def test_CtffindResumePublishesFinishedUnpublishedStepExactlyOnce(
            self,
    ):
        mic = _Mic(
            1,
            "mic_001",
        )

        protocol = (
            _CtffindResumePublishHarness(
                _LogicalMicrographSet(
                    [
                        mic,
                    ],
                    streamClosed=True,
                ),
                existingSteps=[
                    _FinishedCtfStep(
                        "mic_001"
                    ),
                ],
            )
        )

        CistemProtCTFFind.stepsGeneratorStep(
            protocol
        )

        self.assertEqual(
            protocol.scheduledMicNames,
            [],
        )

        self.assertEqual(
            protocol.publishCalls,
            [
                [
                    "mic_001",
                ],
            ],
        )


class _CtffindGeneratorInitialIdsHarness:
    def __init__(self):
        self.finished = False
        self.micDict = None
        self.streamClosed = False
        self.initialStepCalls = 0

    def _defineCtfParamsDict(self):
        pass

    def _insertInitialSteps(self):
        self.initialStepCalls += 1
        return []

    def _restoreProcessedMicsFromPersistentState(self):
        pass

    def _checkNewInput(self):
        if not hasattr(self, "initialIds"):
            raise AssertionError(
                "Generator must initialize initialIds before discovering input."
            )
        self.finished = True

    def _checkNewOutput(self):
        pass

    def _getStreamingSleepOnWait(self):
        return 0

    def _streamingSleepOnWait(self):
        pass


class TestCistemCtffindGeneratorInitialization(unittest.TestCase):
    def test_CtffindGeneratorInitializesInitialIdsBeforeDiscovery(self):
        protocol = _CtffindGeneratorInitialIdsHarness()

        CistemProtCTFFind.stepsGeneratorStep(protocol)

        self.assertEqual(protocol.initialIds, [])
        self.assertEqual(protocol.initialStepCalls, 1)

