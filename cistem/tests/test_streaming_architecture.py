# **************************************************************************
# *
# * Focused regression tests for cisTEM streaming architecture.
# *
# **************************************************************************

import unittest
from collections import OrderedDict

from cistem.protocols.protocol_ctffind import (
    CistemProtCTFFind,
)
from cistem.protocols.protocol_unblur import CistemProtUnblur
from cistem.protocols.protocol_picking import CistemProtFindParticles
from cistem.protocols.protocol_streaming_base import CistemStreamingBase

from .logical_set_fakes import LogicalSetFake


class _Mic:
    def __init__(self, objId, micName):
        self._objId = objId
        self._micName = micName

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return self._micName

    def getFileName(self):
        return self._micName

    def setCTF(self, ctf):
        self._ctf = ctf

    def getCTF(self):
        return getattr(self, '_ctf', None)

    def clone(self):
        return _Mic(
            self._objId,
            self._micName,
        )


class _LogicalMicrographSet(LogicalSetFake):
    pass


class _CtffindStreamingHarness(CistemStreamingBase):
    def __init__(self, inputSet):
        self._inputSet = inputSet
        self.micDict = {}
        self._pendingMics = OrderedDict()
        self._knownMicIds = set()
        self._lastInputId = 0
        self.debugMessages = []

    def getInputMicrographs(self):
        return self._inputSet

    def _getFinishedCtfMicNames(self):
        return CistemProtCTFFind._getFinishedCtfMicNames(self)

    def _getScheduledCtfMicNames(self):
        return CistemProtCTFFind._getScheduledCtfMicNames(self)

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


class _NoJoinStepGuard:
    """Fail loudly if a protocol goes back to the legacy join-step pattern.

    Streaming protocols must not create a WAITING output step that the
    generator unlocks with setStatus(STATUS_NEW); the generator inserts its
    own final step once the stream is exhausted.
    """

    def _getFirstJoinStep(self):
        raise AssertionError(
            "Streaming must not depend on a waiting join step."
        )


class _CtffindInputCheckHarness(
        _NoJoinStepGuard,
        _CtffindStreamingHarness
):
    def __init__(self, inputSet):
        super().__init__(
            inputSet
        )
        self.streamClosed = False
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


class _CtffindOutputCheckHarness(CistemStreamingBase):
    def __init__(self):
        mic = _Mic(
            1,
            "mic_001",
        )

        self.micDict = OrderedDict(
            [(mic.getMicName(), mic)]
        )
        self._pendingMics = OrderedDict()
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

    def getObjId(self):
        return self._mic.getObjId()

    def getMicrograph(self):
        return self._mic


class _LogicalCtfOutput(LogicalSetFake):
    def __init__(self, ctfs):
        super().__init__(ctfs, streamClosed=False)


class _CtffindResumeHarness(CistemStreamingBase):
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

    def _getFinishedCtfMicNames(self):
        return (
            CistemProtCTFFind
            ._getFinishedCtfMicNames(
                self
            )
        )

    def _getScheduledCtfMicNames(self):
        return (
            CistemProtCTFFind
            ._getScheduledCtfMicNames(
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

    def _publishedMicNamesSoFar(self):
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


class _UnblurPointer:
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _LogicalMovie:
    def __init__(self, objId):
        self.objId = objId

    def getObjId(self):
        return self.objId

    def clone(self):
        return _LogicalMovie(self.objId)


class _LogicalMovieSet(LogicalSetFake):
    def __init__(self, movies, streamClosed=True):
        super().__init__(movies, streamClosed=streamClosed)


class _UnblurLogicalInputHarness(CistemStreamingBase):
    def __init__(self, movieSet):
        self.inputMovies = _UnblurPointer(movieSet)
        self.listOfMovies = []
        self.streamClosed = False
        self._lastInputId = 0
        self._knownMovieIds = set()
        self._pendingMovies = {}
        self._discoveredMovieCount = 0

    def debug(self, *args, **kwargs):
        pass


class TestCistemUnblurStreamingArchitecture(unittest.TestCase):
    def test_UnblurStreamingLoadsLogicalMoviesWithoutStorageFilename(self):
        movieSet = _LogicalMovieSet([
            _LogicalMovie(1),
            _LogicalMovie(2),
        ], streamClosed=True)
        protocol = _UnblurLogicalInputHarness(movieSet)

        CistemProtUnblur._loadInputList(protocol)

        self.assertTrue(protocol.streamClosed)
        self.assertEqual([1, 2], [movie.getObjId() for movie in protocol.listOfMovies])
        self.assertGreaterEqual(movieSet.reloads, 1)


class _UnblurInputCheckHarness(_NoJoinStepGuard, CistemStreamingBase):
    def __init__(self, movieSet):
        self.inputMovies = _UnblurPointer(movieSet)
        self.insertedDict = {}
        self.listOfMovies = []
        self.newDeps = []
        self.streamClosed = False
        self.insertedMovieIds = []
        self.updateCalls = 0
        self._lastInputId = 0
        self._knownMovieIds = set()
        self._pendingMovies = {}
        self._discoveredMovieCount = 0

    def debug(self, *args, **kwargs):
        pass

    def _loadInputList(self):
        return CistemProtUnblur._loadInputList(self)

    def _insertNewMoviesSteps(self, insertedDict, inputMovies):
        deps = []
        for movie in inputMovies:
            movieId = movie.getObjId()
            if movieId not in insertedDict:
                self.insertedMovieIds.append(movieId)
                insertedDict[movieId] = movieId
                deps.append(movieId)
        return deps

    def updateSteps(self):
        self.updateCalls += 1


class TestCistemUnblurStreamingInputChecks(unittest.TestCase):
    def test_UnblurStreamingChecksLogicalInputWithoutFilesystemMtime(self):
        movieSet = _LogicalMovieSet([
            _LogicalMovie(1),
            _LogicalMovie(2),
        ], streamClosed=False)
        protocol = _UnblurInputCheckHarness(movieSet)

        CistemProtUnblur._checkNewInput(protocol)

        self.assertEqual([1, 2], protocol.insertedMovieIds)
        self.assertEqual({1: 1, 2: 2}, protocol.insertedDict)
        self.assertFalse(protocol.streamClosed)
        self.assertGreaterEqual(movieSet.reloads, 1)
        self.assertEqual(1, protocol.updateCalls)


class _FinishedMovieStep:
    funcName = "processMovieStep"
    argsStr = '[{"object.id": 1}, false]'

    def isFinished(self):
        return True


class _UnblurOutputCheckHarness(CistemStreamingBase):

    def _getPublishedMovieIds(self):
        return CistemProtUnblur._getPublishedMovieIds(self)

    def _getFinishedMovieIds(self):
        return CistemProtUnblur._getFinishedMovieIds(self)

    def __init__(self):
        movie = _LogicalMovie(1)
        self.listOfMovies = [movie]
        self._pendingMovies = {1: movie}
        self._publishedAnyOutput = False
        self.streamClosed = False
        self.finished = False
        self._steps = [_FinishedMovieStep()]
        self.publishedMovieIds = []

    def _readDoneList(self):
        raise AssertionError(
            "Unblur output publication must not read DONE/all.TXT."
        )

    def _isMovieDone(self, movie):
        raise AssertionError(
            "Unblur output publication must not depend on per-movie DONE sidecars."
        )

    def _updateOutputSets(self, newDone, streamMode):
        self.publishedMovieIds.extend(movie.getObjId() for movie in newDone)

    def debug(self, *args, **kwargs):
        pass


class TestCistemUnblurStreamingCompletion(unittest.TestCase):
    def test_UnblurPublishesFinishedStepWithoutDoneSidecars(self):
        protocol = _UnblurOutputCheckHarness()

        CistemProtUnblur._checkNewOutput(protocol)

        self.assertEqual([1], protocol.publishedMovieIds)


class _UnblurNoDoneProcessingHarness:
    def __init__(self):
        self.convertCIStep = []

    def _getOutputMovieFolder(self, movie):
        return "/tmp"

    def _getMovieDone(self, movie):
        raise AssertionError(
            "Unblur processing must not depend on per-movie DONE sidecars."
        )

    def _filterMovie(self, movie):
        return False


class TestCistemUnblurStreamingProcessing(unittest.TestCase):
    def test_UnblurProcessingDoesNotUsePerMovieDoneSidecars(self):
        from pwem.objects import Movie

        movie = Movie()
        movie.setObjId(1)
        movie.setFileName("/tmp/movie_001.mrc")
        movieDict = movie.getObjDict(includeBasic=True)

        protocol = _UnblurNoDoneProcessingHarness()

        CistemProtUnblur.processMovieStep(protocol, movieDict, False)


class _UnblurCorrectionMovies:
    def __init__(self):
        self.gain = None
        self.dark = None

    def getGain(self):
        return self.gain

    def setGain(self, gain):
        self.gain = gain

    def getDark(self):
        return self.dark

    def setDark(self, dark):
        self.dark = dark


class _UnblurInitialStepHarness:
    def __init__(self):
        self.convertCIStep = []
        self.inputMovies = _UnblurPointer(_UnblurCorrectionMovies())

    def _getExtraPath(self, *parts):
        if parts and parts[0] == "DONE":
            raise AssertionError(
                "Unblur streaming must not create a DONE directory."
            )
        return "/tmp/" + "/".join(parts)

    def _ProtProcessMovies__convertCorrectionImage(self, correctionImage):
        return correctionImage


class TestCistemUnblurStreamingInitialization(unittest.TestCase):
    def test_UnblurInitialConversionDoesNotCreateDoneDirectory(self):
        protocol = _UnblurInitialStepHarness()

        CistemProtUnblur._convertInputStep(protocol)


class _UnblurGeneratorMovies:
    def getSamplingRate(self):
        return 1.5


class _UnblurGeneratorHarness:
    def __init__(self):
        self.inputMovies = _UnblurPointer(_UnblurGeneratorMovies())
        self.insertedFunctions = []
        self.finalDeps = None

    def _insertFunctionStep(self, funcName, *args, **kwargs):
        stepId = len(self.insertedFunctions) + 1
        self.insertedFunctions.append((funcName, args, kwargs, stepId))
        return stepId

    def _insertNewMoviesSteps(self, *args, **kwargs):
        raise AssertionError(
            "Unblur _insertAllSteps must not expand the current movie snapshot."
        )

    def _insertFinalSteps(self, deps):
        self.finalDeps = list(deps)
        return list(deps)


class _UnblurGeneratorLoopHarness(_UnblurGeneratorHarness):
    """Run one generator iteration that already exhausts the input stream."""

    def _restoreProcessedMoviesFromPersistentState(self):
        pass

    def _checkNewInput(self):
        self.newDeps.extend([11, 12])

    def _checkNewOutput(self):
        self.finished = True

    def _getStreamingSleepOnWait(self):
        raise AssertionError("The generator must not sleep once it finished.")


class TestCistemUnblurStreamingGenerator(unittest.TestCase):
    def test_UnblurUsesGeneratorAfterInitialConversion(self):
        protocol = _UnblurGeneratorHarness()

        CistemProtUnblur._insertAllSteps(protocol)

        names = [entry[0] for entry in protocol.insertedFunctions]
        self.assertEqual(["_convertInputStep", "resumableStepGeneratorStep"],
                         names)

        convertCall = next(entry for entry in protocol.insertedFunctions
                           if entry[0] == "_convertInputStep")
        generatorCall = next(entry for entry in protocol.insertedFunctions
                             if entry[0] == "resumableStepGeneratorStep")

        self.assertEqual([convertCall[3]], generatorCall[2].get("prerequisites"))
        # The output step is NOT scheduled upfront with wait=True any more:
        # the generator inserts it once the input stream is exhausted.
        self.assertIsNone(protocol.finalDeps)

    def test_UnblurGeneratorSchedulesOutputStepOnceStreamIsExhausted(self):
        protocol = _UnblurGeneratorLoopHarness()

        CistemProtUnblur.stepsGeneratorStep(protocol)

        names = [entry[0] for entry in protocol.insertedFunctions]
        self.assertEqual(["createOutputStep"], names)

        outputCall = protocol.insertedFunctions[0]
        self.assertEqual([11, 12], outputCall[2].get("prerequisites"))
        self.assertNotIn("wait", outputCall[2])
        self.assertEqual([11, 12], protocol.finalDeps)


class _PendingMovieStep:
    funcName = "processMovieStep"
    argsStr = '[{"object.id": 2}, false]'

    def isFinished(self):
        return False


class _FinishedResumeMovieStep:
    funcName = "processMovieStep"
    argsStr = '[{"object.id": 3}, false]'

    def isFinished(self):
        return True


class _UnblurResumeStateHarness(CistemStreamingBase):
    def __init__(self):
        self.insertedDict = {}
        self._steps = [
            _PendingMovieStep(),
            _FinishedResumeMovieStep(),
        ]
        self.outputMicrographs = _LogicalMovieSet([
            _LogicalMovie(1),
        ], streamClosed=False)
        self.inputMovies = _UnblurPointer(_LogicalMovieSet([
            _LogicalMovie(1),
            _LogicalMovie(2),
            _LogicalMovie(3),
        ], streamClosed=False))

    def _getPublishedMovieIds(self):
        return CistemProtUnblur._getPublishedMovieIds(self)

    def _getScheduledMovieIds(self):
        return CistemProtUnblur._getScheduledMovieIds(self)


class TestCistemUnblurStreamingResumeState(unittest.TestCase):
    def test_UnblurRestoresPublishedAndScheduledMoviesBeforeDiscovery(self):
        protocol = _UnblurResumeStateHarness()

        CistemProtUnblur._restoreProcessedMoviesFromPersistentState(protocol)

        self.assertEqual({1, 2, 3}, set(protocol.insertedDict))


class _ThreadParam:
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _ThreadValidationHarness:
    def __init__(self, threads):
        self.numberOfThreads = _ThreadParam(threads)


class TestCistemUnblurStreamingThreadValidation(unittest.TestCase):
    def test_UnblurGeneratorRejectsOneExecutionWorker(self):
        protocol = _ThreadValidationHarness(2)

        errors = CistemProtUnblur._validateStreamingThreads(protocol)

        self.assertTrue(errors)

    def test_UnblurGeneratorAcceptsTwoExecutionWorkers(self):
        protocol = _ThreadValidationHarness(3)

        errors = CistemProtUnblur._validateStreamingThreads(protocol)

        self.assertEqual([], errors)


class _UnblurFailedSidecarHarness:
    def _getAllFailed(self):
        raise AssertionError(
            "Unblur streaming failure state must not depend on a failed sidecar."
        )


class TestCistemUnblurStreamingFailurePersistence(unittest.TestCase):
    def test_UnblurDoesNotWriteFailedMovieSidecar(self):
        protocol = _UnblurFailedSidecarHarness()

        CistemProtUnblur._writeFailedList(protocol, [_LogicalMovie(1)])


class _UnblurArgFailureHarness:
    def __init__(self):
        self.errors = []

    def getInputMovies(self):
        return object()

    def _createTifLink(self, movie):
        pass

    def _argsUnblur(self, movie):
        raise RuntimeError("argument preparation failed")

    def _getMovieFn(self, movie):
        return movie.getFileName()

    def _getErrorFromUnblurTxt(self, movie, error):
        return str(error)

    def error(self, message):
        self.errors.append(message)


class TestCistemUnblurFailureIsolation(unittest.TestCase):
    def test_UnblurArgBuildingFailureIsToleratedInsteadOfCrashingProtocol(self):
        from pwem.objects import Movie

        movie = Movie()
        movie.setObjId(1)
        movie.setFileName("/tmp/movie_001.mrc")

        protocol = _UnblurArgFailureHarness()

        CistemProtUnblur._processMovie(protocol, movie)

        self.assertEqual(1, len(protocol.errors))
        self.assertIn("argument preparation failed", protocol.errors[0])


class _UnblurMissingShiftsHarness:
    def _getShiftsFn(self, movie):
        return "/tmp/definitely_missing_unblur_shifts.txt"


class TestCistemUnblurMissingShiftsFailure(unittest.TestCase):
    def test_UnblurFailureWithoutShiftsFileIsStillTolerated(self):
        from pwem.objects import Movie

        movie = Movie()
        movie.setObjId(1)
        movie.setFileName("/tmp/movie_001.mrc")
        protocol = _UnblurMissingShiftsHarness()
        originalError = RuntimeError("unblur failed before shifts were written")

        message = CistemProtUnblur._getErrorFromUnblurTxt(
            protocol,
            movie,
            originalError,
        )

        self.assertIs(message, originalError)


class _FindParticlesInputSetGuard:
    def getFileName(self):
        raise AssertionError(
            "FindParticles streaming discovery must not depend on a SQLite/storage filename."
        )


class _FindParticlesInputCheckHarness(_NoJoinStepGuard):
    def __init__(self):
        self.streamClosed = False
        self.micDict = {}
        self.insertedMicNames = []
        self.updateStepsCalls = 0
        self.readyMic = _Mic(1, "mic_001")
        self.inputSet = _FindParticlesInputSetGuard()

    def getInputMicrographs(self):
        return self.inputSet

    def _loadInputList(self):
        return {self.readyMic.getMicName(): self.readyMic}, False

    def _insertNewMicsSteps(self, newMics):
        newMics = list(newMics)
        self.insertedMicNames.extend(mic.getMicName() for mic in newMics)
        return [101 + index for index, _ in enumerate(newMics)]

    def updateSteps(self):
        self.updateStepsCalls += 1


class TestCistemFindParticlesStreamingInputChecks(unittest.TestCase):
    def test_FindParticlesStreamingChecksLogicalInputWithoutFilesystemMtime(self):
        protocol = _FindParticlesInputCheckHarness()

        CistemProtFindParticles._checkNewInput(protocol)

        self.assertEqual(["mic_001"], protocol.insertedMicNames)
        self.assertEqual(1, protocol.updateStepsCalls)
        self.assertFalse(protocol.streamClosed)


class _FinishedPickingStep:
    def __init__(self, funcName, argsStr):
        self.funcName = _StoredValue(funcName)
        self.argsStr = _StoredValue(argsStr)

    def isFinished(self):
        return True


class _PendingPickingStep(_FinishedPickingStep):
    def isFinished(self):
        return False


class _FindParticlesFinishedStepsHarness(CistemStreamingBase):
    def __init__(self):
        self._steps = [
            _FinishedPickingStep(
                "pickMicrographStep",
                '["mic_001", {}]',
            ),
            _FinishedPickingStep(
                "pickMicrographListStep",
                '[["mic_002", "mic_003"], {}]',
            ),
            _PendingPickingStep(
                "pickMicrographStep",
                '["mic_004", {}]',
            ),
        ]


class TestCistemFindParticlesFinishedStepState(unittest.TestCase):
    def test_FindParticlesGetsFinishedMicrographsFromPersistedSteps(self):
        protocol = _FindParticlesFinishedStepsHarness()

        finished = CistemProtFindParticles._getFinishedPickingMicNames(protocol)

        self.assertEqual(
            {"mic_001", "mic_002", "mic_003"},
            finished,
        )


class _FindParticlesOutputCheckHarness(CistemStreamingBase):
    def __init__(self):
        mic = _Mic(1, "mic_001")
        self.micDict = OrderedDict([(mic.getMicName(), mic)])
        self._pendingMics = OrderedDict()
        self._micsWithoutCtf = OrderedDict()
        self.streamClosed = False
        self.finished = False
        self._steps = [
            _FinishedPickingStep(
                "pickMicrographStep",
                '["mic_001", {}]',
            ),
        ]
        self.publishedMicNames = set()
        self.publishCalls = []
        self.sleepCalls = 0

    def _readDoneList(self):
        raise AssertionError(
            "FindParticles output publication must not depend on DONE/all.TXT."
        )

    def _isMicDone(self, mic):
        raise AssertionError(
            "FindParticles output publication must not depend on DONE/mic_*.TXT."
        )

    def _getFinishedPickingMicNames(self):
        return CistemProtFindParticles._getFinishedPickingMicNames(self)

    def _updateOutputCoordSet(self, micList, streamMode):
        micList = list(micList)
        self.publishCalls.append([mic.getMicName() for mic in micList])
        self.publishedMicNames.update(mic.getMicName() for mic in micList)
        return micList

    def _streamingSleepOnWait(self):
        self.sleepCalls += 1

    def debug(self, *args, **kwargs):
        pass


class TestCistemFindParticlesStreamingCompletion(unittest.TestCase):
    def test_FindParticlesPublishesFinishedStepWithoutDoneSidecars(self):
        protocol = _FindParticlesOutputCheckHarness()

        CistemProtFindParticles._checkNewOutput(protocol)
        CistemProtFindParticles._checkNewOutput(protocol)

        self.assertEqual(
            [["mic_001"]],
            protocol.publishCalls,
        )
        self.assertEqual(
            {"mic_001"},
            protocol.publishedMicNames,
        )


class _FindParticlesEmptyCoordHarness(CistemStreamingBase):
    """A micrograph that was picked but yielded no coordinate at all."""

    def __init__(self):
        self.micDict = OrderedDict()
        self._pendingMics = OrderedDict()
        self._micsWithoutCtf = OrderedDict()
        self._ctfByMicName = {}
        self._knownMicIds = set()
        self._knownCtfIds = set()
        self._lastMicId = 0
        self._lastCtfId = 0
        self.streamClosed = False
        self.scheduledMicNames = []
        self._inputMics = _FindParticlesResumeLogicalSet([
            _Mic(1, "mic_001"),
        ])
        self.ctfRelations = _UnblurPointer(_FindParticlesResumeLogicalSet([
            _LogicalCtf(_Mic(1, "mic_001")),
        ]))
        # It was picked, so its step is in the graph, but it left no row
        # behind in the output coordinates.
        self._steps = [
            _FinishedPickingStep("pickMicrographStep", '["mic_001", {}]'),
        ]
        self.outputCoordinates = _LogicalCoordinateOutput([])

    def getInputMicrographs(self):
        return self._inputMics

    def _getPublishedPickingMicIds(self):
        return CistemProtFindParticles._getPublishedPickingMicIds(self)

    def _getScheduledPickingMicNames(self):
        return CistemProtFindParticles._getScheduledPickingMicNames(self)

    def _loadInputList(self):
        return CistemProtFindParticles._loadInputList(self)

    def _loadMics(self, micSet):
        return CistemProtFindParticles._loadMics(self, micSet)

    def _loadCTFs(self, ctfSet):
        return CistemProtFindParticles._loadCTFs(self, ctfSet)

    def _loadSet(self, inputSet, SetClass, getKeyFunc, watermarkAttr,
                 knownIds):
        return CistemProtFindParticles._loadSet(
            self, inputSet, SetClass, getKeyFunc, watermarkAttr, knownIds)

    def _insertNewMicsSteps(self, newMics):
        newMics = list(newMics)
        self.scheduledMicNames.extend(mic.getMicName() for mic in newMics)

        for mic in newMics:
            self.micDict[mic.getMicName()] = mic

        return []

    def updateSteps(self):
        pass

    def debug(self, *args, **kwargs):
        pass


class TestCistemFindParticlesPublishedStatePersistence(unittest.TestCase):
    def test_FindParticlesDoesNotRepickMicrographThatYieldedNoCoordinates(self):
        # The output Set cannot answer this one - a micrograph with zero
        # coordinates leaves no row - so the step graph has to, which is
        # why no separate published-name list has to be persisted.
        protocol = _FindParticlesEmptyCoordHarness()

        CistemProtFindParticles._restoreProcessedMicsFromPersistentState(
            protocol)
        CistemProtFindParticles._checkNewInput(protocol)

        self.assertEqual([], protocol.scheduledMicNames)
        self.assertEqual({"mic_001"}, set(protocol.micDict))


class _FindParticlesNoDoneProcessingHarness:
    def __init__(self):
        self.micDict = {
            "mic_001": _Mic(1, "mic_001"),
            "mic_002": _Mic(2, "mic_002"),
        }
        self.pickedMicNames = []
        self.batchPickedMicNames = []

    def isContinued(self):
        return True

    def _getMicDone(self, mic):
        raise AssertionError(
            "FindParticles processing must not depend on per-micrograph DONE sidecars."
        )

    def _pickMicrograph(self, mic, *args):
        self.pickedMicNames.append(mic.getMicName())

    def _pickMicrographList(self, micList, *args):
        self.batchPickedMicNames.extend(mic.getMicName() for mic in micList)

    def info(self, *args, **kwargs):
        pass


class TestCistemFindParticlesStreamingProcessing(unittest.TestCase):
    def test_FindParticlesSingleProcessingDoesNotUseDoneSidecar(self):
        protocol = _FindParticlesNoDoneProcessingHarness()

        CistemProtFindParticles.pickMicrographStep(
            protocol,
            "mic_001",
        )

        self.assertEqual(
            ["mic_001"],
            protocol.pickedMicNames,
        )

    def test_FindParticlesBatchProcessingDoesNotUseDoneSidecars(self):
        protocol = _FindParticlesNoDoneProcessingHarness()

        CistemProtFindParticles.pickMicrographListStep(
            protocol,
            ["mic_001", "mic_002"],
        )

        self.assertEqual(
            ["mic_001", "mic_002"],
            protocol.batchPickedMicNames,
        )


class _FindParticlesFailedSidecarHarness:
    def _getAllFailed(self):
        raise AssertionError(
            "FindParticles failure state must not depend on FAILED_all.TXT."
        )


class TestCistemFindParticlesFailurePersistence(unittest.TestCase):
    def test_FindParticlesDoesNotWriteFailedMicrographSidecar(self):
        protocol = _FindParticlesFailedSidecarHarness()

        CistemProtFindParticles._writeFailedList(
            protocol,
            [_Mic(1, "mic_001")],
        )


class _FindParticlesStreamingInput:
    def isStreamOpen(self):
        return True


class _FindParticlesStreamingInsertHarness:
    streamingBatchSize = 1

    def __init__(self):
        self.inputStreaming = False
        self.insertCalls = []
        self.legacyInitialCalls = 0
        self.legacyLoadCalls = 0
        self.legacyFinalCalls = 0

    def getInputMicrographs(self):
        return _FindParticlesStreamingInput()

    def resumableStepGeneratorStep(self, timestamp):
        pass

    def _insertFunctionStep(self, funcName, *args, **kwargs):
        if callable(funcName):
            funcName = funcName.__name__

        self.insertCalls.append({
            "funcName": funcName,
            "args": args,
            "kwargs": kwargs,
        })
        return len(self.insertCalls)

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

    def _getAllDone(self):
        raise AssertionError(
            "FindParticles streaming orchestration must not create DONE/all.TXT."
        )


class TestCistemFindParticlesGeneratorOrchestration(unittest.TestCase):
    def test_FindParticlesStreamingUsesSingleResumableGenerator(self):
        protocol = _FindParticlesStreamingInsertHarness()

        CistemProtFindParticles._insertAllSteps(protocol)

        self.assertTrue(protocol.inputStreaming)
        self.assertEqual(0, protocol.legacyInitialCalls)
        self.assertEqual(0, protocol.legacyLoadCalls)
        self.assertEqual(0, protocol.legacyFinalCalls)

        self.assertEqual(1, len(protocol.insertCalls))
        self.assertEqual(
            "resumableStepGeneratorStep",
            protocol.insertCalls[0]["funcName"],
        )
        self.assertFalse(
            protocol.insertCalls[0]["kwargs"].get("needsGPU", True)
        )


class _ScheduledPickingStep:
    def __init__(self, funcName, argsStr, finished=False):
        self.funcName = _StoredValue(funcName)
        self.argsStr = _StoredValue(argsStr)
        self._finished = finished

    def isFinished(self):
        return self._finished


class _FindParticlesResumeLogicalSet(LogicalSetFake):
    def __init__(self, items, streamClosed=True):
        super().__init__(items, streamClosed=streamClosed)


class _LogicalCtf:
    def __init__(self, mic):
        self._mic = mic

    def getObjId(self):
        return self._mic.getObjId()

    def getMicrograph(self):
        return self._mic

    def clone(self):
        return _LogicalCtf(self._mic)


class _PublishedCoordinate:
    """A coordinate only tells which micrograph it came from."""

    def __init__(self, micId):
        self._micId = micId

    def getObjId(self):
        return self._micId

    def getMicId(self):
        return self._micId


class _LogicalCoordinateOutput(LogicalSetFake):
    def __init__(self, micIds):
        super().__init__([_PublishedCoordinate(micId) for micId in micIds],
                         streamClosed=False)

    def getUniqueValues(self, attributes, where=None):
        if attributes == '_micId':
            return [item.getMicId() for item in self._items]

        return super().getUniqueValues(attributes, where=where)


class _FindParticlesResumeHarness(CistemStreamingBase):
    def __init__(self):
        self.micDict = OrderedDict()
        self._pendingMics = OrderedDict()
        self._micsWithoutCtf = OrderedDict()
        self._ctfByMicName = {}
        self._knownMicIds = set()
        self._knownCtfIds = set()
        self._lastMicId = 0
        self._lastCtfId = 0
        self.streamClosed = False
        self.finished = False
        self.initialIds = []
        self.scheduledMicNames = []
        self._inputMics = _FindParticlesResumeLogicalSet([
            _Mic(1, "mic_001"),
            _Mic(2, "mic_002"),
            _Mic(3, "mic_003"),
        ])
        self._steps = [
            _ScheduledPickingStep(
                "pickMicrographStep",
                '["mic_002", {}]',
                finished=False,
            ),
        ]
        # mic_001 was published by a previous run, mic_002 only scheduled.
        self.outputCoordinates = _LogicalCoordinateOutput([1])
        # Every micrograph already has its CTF in this scenario.
        self.ctfRelations = _UnblurPointer(_FindParticlesResumeLogicalSet([
            _LogicalCtf(_Mic(1, "mic_001")),
            _LogicalCtf(_Mic(2, "mic_002")),
            _LogicalCtf(_Mic(3, "mic_003")),
        ]))

    def getInputMicrographs(self):
        return self._inputMics

    def _getPublishedPickingMicIds(self):
        return CistemProtFindParticles._getPublishedPickingMicIds(self)

    def _getScheduledPickingMicNames(self):
        return CistemProtFindParticles._getScheduledPickingMicNames(self)

    def _restoreProcessedMicsFromPersistentState(self):
        return CistemProtFindParticles._restoreProcessedMicsFromPersistentState(self)

    def _loadInputList(self):
        return CistemProtFindParticles._loadInputList(self)

    def _loadMics(self, micSet):
        return CistemProtFindParticles._loadMics(self, micSet)

    def _loadCTFs(self, ctfSet):
        return CistemProtFindParticles._loadCTFs(self, ctfSet)

    def _loadSet(self, inputSet, SetClass, getKeyFunc, watermarkAttr,
                 knownIds):
        return CistemProtFindParticles._loadSet(
            self, inputSet, SetClass, getKeyFunc, watermarkAttr, knownIds)

    def _insertNewMicsSteps(self, newMics):
        newMics = list(newMics)
        self.scheduledMicNames.extend(mic.getMicName() for mic in newMics)
        for mic in newMics:
            self.micDict[mic.getMicName()] = mic
        return []

    def updateSteps(self):
        pass

    def debug(self, *args, **kwargs):
        pass


class TestCistemFindParticlesGeneratorResumeSafety(unittest.TestCase):
    def test_FindParticlesResumeDoesNotReschedulePublishedOrScheduledMicrographs(self):
        protocol = _FindParticlesResumeHarness()

        CistemProtFindParticles._restoreProcessedMicsFromPersistentState(protocol)
        CistemProtFindParticles._checkNewInput(protocol)

        self.assertEqual(
            ["mic_003"],
            protocol.scheduledMicNames,
        )
        # micDict now holds only what still has to be published, so the
        # micrograph a previous run already published stays out of it.
        self.assertEqual(
            {"mic_002", "mic_003"},
            set(protocol.micDict),
        )


class _FindParticlesThreadParam:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _FindParticlesThreadInput:
    def __init__(self, streamOpen):
        self._streamOpen = streamOpen

    def isStreamOpen(self):
        return self._streamOpen


class _FindParticlesThreadValidationHarness:
    def __init__(self, threads, streamOpen):
        self.numberOfThreads = _FindParticlesThreadParam(threads)
        self._inputMics = _FindParticlesThreadInput(streamOpen)

    def getInputMicrographs(self):
        return self._inputMics


class TestCistemFindParticlesStreamingThreadValidation(unittest.TestCase):
    def test_FindParticlesStreamingRejectsOneExecutionWorker(self):
        protocol = _FindParticlesThreadValidationHarness(
            threads=2,
            streamOpen=True,
        )

        errors = CistemProtFindParticles._validateStreamingThreads(protocol)

        self.assertTrue(errors)

    def test_FindParticlesStreamingAcceptsTwoExecutionWorkers(self):
        protocol = _FindParticlesThreadValidationHarness(
            threads=3,
            streamOpen=True,
        )

        errors = CistemProtFindParticles._validateStreamingThreads(protocol)

        self.assertEqual([], errors)

    def test_FindParticlesNonStreamingStillAcceptsOneThread(self):
        protocol = _FindParticlesThreadValidationHarness(
            threads=1,
            streamOpen=False,
        )

        errors = CistemProtFindParticles._validateStreamingThreads(protocol)

        self.assertEqual([], errors)



class _ResumedStepsHarness(CistemStreamingBase):
    """Work finished before a Continue lives in _prevSteps, not _steps."""

    def __init__(self, prevSteps, steps=None):
        self._prevSteps = prevSteps
        self._steps = steps or []


class TestCistemResumedStepGraphVisibility(unittest.TestCase):
    def test_FinishedStepsRestoredOnResumeAreNotRescheduled(self):
        # Reading only _steps would hide everything a previous execution
        # already finished, so Continue would redo that work.
        protocol = _ResumedStepsHarness(
            prevSteps=[
                _FinishedPickingStep(
                    "pickMicrographListStep",
                    '[["mic_001", "mic_002"], {}]',
                ),
            ],
            steps=[
                _FinishedPickingStep(
                    "pickMicrographStep",
                    '["mic_003", {}]',
                ),
            ],
        )

        finished = CistemProtFindParticles._getFinishedPickingMicNames(protocol)

        self.assertEqual({"mic_001", "mic_002", "mic_003"}, finished)

    def test_ResumedStepsAreCountedOnlyOnceWhenSharedBetweenLists(self):
        sharedStep = _FinishedPickingStep(
            "pickMicrographStep",
            '["mic_001", {}]',
        )

        protocol = _ResumedStepsHarness(
            prevSteps=[sharedStep],
            steps=[sharedStep],
        )

        self.assertEqual(
            1,
            len(list(protocol._iterKnownStreamingSteps())),
        )

    def test_UnblurResumeAlsoSeesStepsRestoredFromPreviousRun(self):
        protocol = _ResumedStepsHarness(
            prevSteps=[_FinishedMovieStep()],
        )

        self.assertEqual(
            {1},
            CistemProtUnblur._getFinishedMovieIds(protocol),
        )


class TestCistemCtffindStreamingThreadValidation(unittest.TestCase):
    def test_CtffindGeneratorRejectsOneExecutionWorker(self):
        protocol = _ThreadValidationHarness(2)

        errors = CistemProtCTFFind._validateStreamingThreads(protocol)

        self.assertTrue(errors)

    def test_CtffindGeneratorAcceptsTwoExecutionWorkers(self):
        protocol = _ThreadValidationHarness(3)

        errors = CistemProtCTFFind._validateStreamingThreads(protocol)

        self.assertEqual([], errors)


class _CostTrackingHarness(_CtffindStreamingHarness):
    """Ctffind discovery over a Set that reports what it hydrated."""

    def __init__(self, inputSet):
        super().__init__(inputSet)
        self._steps = []
        self.insertedMicNames = []
        self.updateStepsCalls = 0

    def _getFinishedCtfMicNames(self):
        return CistemProtCTFFind._getFinishedCtfMicNames(self)

    def _getScheduledCtfMicNames(self):
        return CistemProtCTFFind._getScheduledCtfMicNames(self)

    def _loadInputList(self):
        return CistemProtCTFFind._loadSet(
            self, self._inputSet, None, lambda mic: mic.getMicName())

    def _insertNewMicsSteps(self, newMics):
        newMics = list(newMics)
        self.insertedMicNames.extend(mic.getMicName() for mic in newMics)

        for mic in newMics:
            self.micDict[mic.getMicName()] = mic

        return []

    def updateSteps(self):
        self.updateStepsCalls += 1


class TestCistemStreamingDiscoveryCost(unittest.TestCase):
    """A poll must cost what just arrived, not everything seen so far."""

    @staticmethod
    def _mics(firstId, count):
        return [_Mic(micId, "mic_%03d" % micId)
                for micId in range(firstId, firstId + count)]

    def test_CtffindPollOnlyHydratesTheMicrographsThatJustArrived(self):
        inputSet = _LogicalMicrographSet(self._mics(1, 500))
        protocol = _CostTrackingHarness(inputSet)

        CistemProtCTFFind._checkNewInput(protocol)

        self.assertEqual(500, inputSet.hydratedItems)

        inputSet.addItems(self._mics(501, 3))
        hydratedBefore = inputSet.hydratedItems

        CistemProtCTFFind._checkNewInput(protocol)

        self.assertEqual(3, inputSet.hydratedItems - hydratedBefore)
        self.assertEqual(0, inputSet.fullScans)

    def test_CtffindIdlePollHydratesNothingAtAll(self):
        inputSet = _LogicalMicrographSet(self._mics(1, 500))
        protocol = _CostTrackingHarness(inputSet)

        CistemProtCTFFind._checkNewInput(protocol)
        hydratedBefore = inputSet.hydratedItems

        CistemProtCTFFind._checkNewInput(protocol)

        self.assertEqual(hydratedBefore, inputSet.hydratedItems)

    def test_UnblurPollOnlyHydratesTheMoviesThatJustArrived(self):
        movies = [_LogicalMovie(movieId) for movieId in range(1, 501)]
        movieSet = _LogicalMovieSet(movies, streamClosed=False)
        protocol = _UnblurInputCheckHarness(movieSet)

        CistemProtUnblur._checkNewInput(protocol)

        self.assertEqual(500, movieSet.hydratedItems)

        movieSet.addItems([_LogicalMovie(501), _LogicalMovie(502)])
        hydratedBefore = movieSet.hydratedItems

        CistemProtUnblur._checkNewInput(protocol)

        self.assertEqual(2, movieSet.hydratedItems - hydratedBefore)
        self.assertEqual(0, movieSet.fullScans)
        self.assertEqual([501, 502], protocol.insertedMovieIds[-2:])

    def test_StepGraphScanDoesNotReparseStepsItAlreadyRead(self):
        protocol = _CostTrackingHarness(_LogicalMicrographSet([]))

        parsed = []
        original = CistemStreamingBase._parseStepArgKeys

        def countingParse(step, dictField, keyType):
            parsed.append(step)
            return original(step, dictField, keyType)

        # A plain function, not staticmethod(): an instance attribute is
        # never bound, and a staticmethod object is only callable itself
        # from Python 3.10 on.
        protocol._parseStepArgKeys = countingParse

        protocol._steps = [
            _FinishedCtfStep("mic_%03d" % micId) for micId in range(1, 201)
        ]

        self.assertEqual(200, len(protocol._getFinishedCtfMicNames()))
        self.assertEqual(200, len(parsed))

        protocol._steps.append(_FinishedCtfStep("mic_201"))

        self.assertEqual(201, len(protocol._getFinishedCtfMicNames()))
        # Only the new step is parsed again, not the 200 already read.
        self.assertEqual(201, len(parsed))

    def test_ResumeStartsAboveWhatIsAlreadyPublishedButKeepsTheGaps(self):
        inputSet = _LogicalMicrographSet(self._mics(1, 100))
        protocol = _CostTrackingHarness(inputSet)

        # A previous run published everything except mic_042.
        publishedIds = set(range(1, 101)) - {42}

        watermark, gapIds = protocol._resumeWatermarkWithGaps(
            inputSet, publishedIds)

        self.assertEqual(100, watermark)
        self.assertEqual({42}, gapIds)
        self.assertEqual(0, inputSet.hydratedItems)

    def test_StepGraphScanSeesAPendingStepFinishAfterTheListIsRebuilt(self):
        protocol = _CostTrackingHarness(_LogicalMicrographSet([]))
        protocol._steps = [
            _ScheduledPickingStep(
                "estimateCtfStep", '["mic_001", {}]', finished=False),
        ]

        self.assertEqual(set(), protocol._getFinishedCtfMicNames())

        # pyworkflow can rebuild the step list from the database, handing
        # back new objects for the same positions.
        protocol._steps = [
            _ScheduledPickingStep(
                "estimateCtfStep", '["mic_001", {}]', finished=True),
        ]

        self.assertEqual({"mic_001"}, protocol._getFinishedCtfMicNames())
