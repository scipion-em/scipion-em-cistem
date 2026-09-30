# **************************************************************************
# *
# * Regression tests for CISTEM streaming/Continue behaviour.
# *
# **************************************************************************

import os
import tempfile
import unittest
from unittest.mock import patch

from cistem.protocols.protocol_picking import CistemProtFindParticles
from cistem.protocols.protocol_unblur import CistemProtUnblur
from cistem.protocols.protocol_ctffind import CistemProtCTFFind


class _Mic:
    def __init__(self, obj_id=1, ctf=None):
        self._obj_id = obj_id
        self._ctf = ctf

    def getObjId(self):
        return self._obj_id

    def getFileName(self):
        return "mic_001.mrc"

    def getMicName(self):
        return "mic_001"

    def getCTF(self):
        return self._ctf


class _Ctf:
    def getPhaseShift(self):
        return 0.0

    def getDefocusU(self):
        return 10000.0

    def getDefocusV(self):
        return 11000.0

    def getDefocusAngle(self):
        return 15.0


class _PickingHarness:
    def __init__(self, tmp_dir):
        self.tmp_dir = tmp_dir
        self.ctfDict = {"mic_001": _Ctf()}
        self.pickType = 0
        self.errors = []
        self.failed_ids = []

    def _getTmpPath(self, *parts):
        return os.path.join(self.tmp_dir, *parts)

    def _getLogFn(self, mic):
        return os.path.join(self.tmp_dir, "missing_picker.log")

    def _getStackFn(self, mic):
        return os.path.join(self.tmp_dir, "picker.mrc")

    def _getPltFn(self, mic):
        return os.path.join(self.tmp_dir, "picker.plt")

    def _getProgram(self):
        return "find_particles"

    def _getArgsStr(self):
        return "%(micName)s"

    def runJob(self, *args, **kwargs):
        raise RuntimeError("simulated find_particles failure")

    def _getErrorFromPickerTxt(self, mic, error):
        return CistemProtFindParticles._getErrorFromPickerTxt(
            self, mic, error
        )

    def _convertMic(self, mic):
        return CistemProtFindParticles._convertMic(self, mic)

    def _getMicCtf(self, mic):
        return CistemProtFindParticles._getMicCtf(self, mic)

    def _writeFailedList(self, mics):
        self.failed_ids.extend(
            mic.getObjId() for mic in mics
        )

    def error(self, message):
        self.errors.append(message)


class TestCistemStreamingRegression(unittest.TestCase):
    def testPickingFailureWithoutLogIsStillTolerated(self):
        with tempfile.TemporaryDirectory() as tmp:
            protocol = _PickingHarness(tmp)
            mic = _Mic()

            CistemProtFindParticles._pickMicrographStep(
                protocol,
                [mic],
                {},
            )

            self.assertEqual(1, len(protocol.errors))
            self.assertIn(
                "simulated find_particles failure",
                protocol.errors[0],
            )


    def testCtffindFailureWithoutLogIsStillTolerated(self):
        with tempfile.TemporaryDirectory() as tmp:
            mic_file = os.path.join(tmp, "mic_001.mrc")
            with open(mic_file, "w") as handle:
                handle.write("test")

            class _CtfMic:
                def getObjId(self):
                    return 1

                def getFileName(self):
                    return mic_file

            class _Program:
                def getCommand(self, **kwargs):
                    return "ctffind", "args"

            class _CtffindHarness:
                usePowerSpectra = False

                def __init__(self):
                    self._ctfProgram = _Program()
                    self.errors = []

                def _getTmpPath(self, *parts):
                    return os.path.join(tmp, *parts)

                def _getCtfOutPath(self, mic):
                    return os.path.join(tmp, "missing_ctf.txt")

                def _getPsdPath(self, mic):
                    return os.path.join(tmp, "ctf.mrc")

                def runJob(self, *args, **kwargs):
                    raise RuntimeError("simulated ctffind failure")

                def _getErrorFromCtffindTxt(self, mic, error):
                    return CistemProtCTFFind._getErrorFromCtffindTxt(
                        self, mic, error
                    )

                def error(self, message):
                    self.errors.append(message)

            protocol = _CtffindHarness()

            CistemProtCTFFind._doCtfEstimation(
                protocol,
                _CtfMic(),
            )

            self.assertEqual(1, len(protocol.errors))
            self.assertIn(
                "simulated ctffind failure",
                protocol.errors[0],
            )

    def testCtffindMissingMicrographIsLoggedInsteadOfCrashingProtocol(self):
        # Regression test: the missing-micrograph existence check (and
        # the mrc conversion right after it) used to sit OUTSIDE the
        # try/except, raising uncaught. For streamingBatchSize == 1
        # (ctffind's default), _estimateCTF is called directly from
        # estimateCtfStep with no exception boundary of its own around
        # it, so this would crash the whole protocol for a single
        # missing/corrupted micrograph, instead of being logged and
        # skipped like every other failure in this function.
        with tempfile.TemporaryDirectory() as tmp:
            missing_mic_file = os.path.join(tmp, "does_not_exist.mrc")

            class _CtfMic:
                def getObjId(self):
                    return 1

                def getFileName(self):
                    return missing_mic_file

            class _Program:
                def getCommand(self, **kwargs):
                    raise AssertionError(
                        "The ctffind command must not be built for a "
                        "micrograph whose input file does not exist."
                    )

            class _CtffindHarness:
                usePowerSpectra = False

                def __init__(self):
                    self._ctfProgram = _Program()
                    self.errors = []

                def _getTmpPath(self, *parts):
                    return os.path.join(tmp, *parts)

                def _getCtfOutPath(self, mic):
                    return os.path.join(tmp, "missing_ctf.txt")

                def _getPsdPath(self, mic):
                    return os.path.join(tmp, "ctf.mrc")

                def runJob(self, *args, **kwargs):
                    raise AssertionError(
                        "runJob must not be called for a missing "
                        "micrograph."
                    )

                def _getErrorFromCtffindTxt(self, mic, error):
                    return CistemProtCTFFind._getErrorFromCtffindTxt(
                        self, mic, error
                    )

                def error(self, message):
                    self.errors.append(message)

            protocol = _CtffindHarness()

            CistemProtCTFFind._doCtfEstimation(
                protocol,
                _CtfMic(),
            )

            self.assertEqual(1, len(protocol.errors))
            self.assertIn(
                "Missing input micrograph",
                protocol.errors[0],
            )


    def testUnblurFailureWithoutShiftsFileIsStillTolerated(self):
        with tempfile.TemporaryDirectory() as tmp:
            class _InputMovies:
                def getSamplingRate(self):
                    return 1.0

            class _Movie:
                pass

            class _UnblurHarness:
                def __init__(self):
                    self.errors = []
                    self._args = ""

                def getInputMovies(self):
                    return _InputMovies()

                def _createTifLink(self, movie):
                    pass

                def _argsUnblur(self, movie):
                    self._args = "args"

                def _getProgram(self):
                    return "unblur"

                def runJob(self, *args, **kwargs):
                    raise RuntimeError("simulated unblur failure")

                def _getMovieFn(self, movie):
                    return os.path.join(tmp, "movie_001.mrc")

                def _getShiftsFn(self, movie):
                    return os.path.join(tmp, "missing_shifts.txt")

                def _getErrorFromUnblurTxt(self, movie, error):
                    return CistemProtUnblur._getErrorFromUnblurTxt(
                        self, movie, error
                    )

                def error(self, message):
                    self.errors.append(message)

            protocol = _UnblurHarness()

            with patch(
                "cistem.protocols.protocol_unblur.Plugin.getEnviron",
                return_value={},
            ):
                CistemProtUnblur._processMovie(
                    protocol,
                    _Movie(),
                )

            self.assertEqual(1, len(protocol.errors))
            self.assertIn(
                "simulated unblur failure",
                protocol.errors[0],
            )

    def testUnblurArgBuildingFailureIsToleratedInsteadOfCrashingProtocol(self):
        # Regression test: _createTifLink/_argsUnblur used to run BEFORE
        # the try/except. processMovieStep (the pwem base class step
        # calling _processMovie) has no exception boundary of its own,
        # so a single movie failing here (e.g. a bad tiff link, or a
        # missing acquisition attribute while building the unblur
        # arguments) would crash the whole protocol instead of being
        # logged and skipped like every other failure in this function.
        with tempfile.TemporaryDirectory() as tmp:
            class _InputMovies:
                def getSamplingRate(self):
                    return 1.0

            class _Movie:
                pass

            class _UnblurHarness:
                def __init__(self):
                    self.errors = []
                    self._args = ""

                def getInputMovies(self):
                    return _InputMovies()

                def _createTifLink(self, movie):
                    pass

                def _argsUnblur(self, movie):
                    raise KeyError("acquisition")

                def _getProgram(self):
                    raise AssertionError(
                        "The unblur program must not be resolved when "
                        "argument building already failed."
                    )

                def runJob(self, *args, **kwargs):
                    raise AssertionError(
                        "runJob must not be called when argument "
                        "building already failed."
                    )

                def _getMovieFn(self, movie):
                    return os.path.join(tmp, "movie_001.mrc")

                def _getShiftsFn(self, movie):
                    return os.path.join(tmp, "missing_shifts.txt")

                def _getErrorFromUnblurTxt(self, movie, error):
                    return CistemProtUnblur._getErrorFromUnblurTxt(
                        self, movie, error
                    )

                def error(self, message):
                    self.errors.append(message)

            protocol = _UnblurHarness()

            CistemProtUnblur._processMovie(
                protocol,
                _Movie(),
            )

            self.assertEqual(1, len(protocol.errors))


    def testPickingFailureIsRecordedBeforeContinuing(self):
        with tempfile.TemporaryDirectory() as tmp:
            protocol = _PickingHarness(tmp)
            protocol.failed_ids = []

            def _record_failed(mics):
                protocol.failed_ids.extend(
                    mic.getObjId() for mic in mics
                )

            protocol._writeFailedList = _record_failed
            mic = _Mic()

            CistemProtFindParticles._pickMicrographStep(
                protocol,
                [mic],
                {},
            )

            self.assertEqual(
                [1],
                protocol.failed_ids,
                "A tolerated picking failure must remain distinguishable "
                "from a valid micrograph with zero picked particles.",
            )


    def testPickMicrographStepConvertsMicLazilyWhenConvertInputStepMissedIt(self):
        # Regression test: convertInputStep only runs ONCE, near the
        # start of the protocol, as an initial step - but a streaming
        # input Set keeps growing afterwards. A micrograph that arrives
        # later was never converted to mrc by that one-shot snapshot,
        # so _pickMicrographStep must convert it lazily itself instead
        # of silently failing against a missing input file for every
        # micrograph that streams in after that snapshot.
        with tempfile.TemporaryDirectory() as tmp:
            protocol = _PickingHarness(tmp)
            mic = _Mic()

            outMic = os.path.join(tmp, "mic_001.mrc")
            self.assertFalse(os.path.lexists(outMic))

            CistemProtFindParticles._pickMicrographStep(
                protocol,
                [mic],
                {},
            )

            # lexists (not exists) because the fake mic's source file
            # does not really exist on disk here - only that the
            # conversion (a symlink, for a .mrc source) was attempted
            # is being asserted, not that its target is readable.
            self.assertTrue(
                os.path.lexists(outMic),
                "_pickMicrographStep must convert a micrograph lazily "
                "when convertInputStep never processed it.",
            )

    def testGetMicCtfPrefersFreshMicCtfOverStaleCtfDictSnapshot(self):
        # Regression test: streaming input already has its CTF attached
        # via mic.setCTF() (kept fresh by _loadInputList on every
        # poll). self.ctfDict is a one-shot snapshot built once in
        # convertInputStep and never refreshed - a micrograph whose CTF
        # was computed after that snapshot would be missing from (or
        # stale in) ctfDict, so the fresh mic.getCTF() must win.
        freshCtf = _Ctf()
        staleCtf = _Ctf()
        mic = _Mic(ctf=freshCtf)

        protocol = _PickingHarness("unused")
        protocol.ctfDict = {"mic_001": staleCtf}

        ctf = CistemProtFindParticles._getMicCtf(protocol, mic)

        self.assertIs(freshCtf, ctf)

    def testGetMicCtfFallsBackToSnapshotForNonStreamingMics(self):
        # Non-streaming micrographs never go through _loadInputList, so
        # they never have a CTF attached via setCTF() - convertInputStep's
        # one-shot snapshot remains the only source for them.
        mic = _Mic(ctf=None)
        snapshotCtf = _Ctf()

        protocol = _PickingHarness("unused")
        protocol.ctfDict = {"mic_001": snapshotCtf}

        ctf = CistemProtFindParticles._getMicCtf(protocol, mic)

        self.assertIs(snapshotCtf, ctf)

    def testPickMicrographStepSkipsMicWithNoCtfAvailableAnywhereInsteadOfCrashing(self):
        # Regression test: the original self.ctfDict[mic.getMicName()]
        # was a direct, unguarded dict index - a KeyError for a
        # micrograph missing from both the fresh mic.getCTF() and the
        # stale snapshot would propagate uncaught out of
        # _pickMicrographStep, crashing the whole batch instead of
        # being isolated per-micrograph like every other picking
        # failure.
        with tempfile.TemporaryDirectory() as tmp:
            protocol = _PickingHarness(tmp)
            protocol.ctfDict = {}
            mic = _Mic(ctf=None)

            CistemProtFindParticles._pickMicrographStep(
                protocol,
                [mic],
                {},
            )

            self.assertEqual(1, len(protocol.errors))
            self.assertEqual([1], protocol.failed_ids)

    def testPickingStreamingLoadsLogicalSetsWithoutStorageFilename(self):
        from collections import OrderedDict

        class _Pointer:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        class _LogicalMic:
            def __init__(self, objId, name):
                self.objId = objId
                self.name = name
                self.ctf = None

            def getObjId(self):
                return self.objId

            def getMicName(self):
                return self.name

            def setCTF(self, ctf):
                self.ctf = ctf

            def clone(self):
                clone = _LogicalMic(self.objId, self.name)
                clone.ctf = self.ctf
                return clone

        class _LogicalCtf:
            def __init__(self, mic):
                self.mic = mic

            def getMicrograph(self):
                return self.mic

            def clone(self):
                return _LogicalCtf(self.mic)

        class _LogicalSet:
            def __init__(self, items, closed=True):
                self.items = list(items)
                self.closed = closed
                self.reloads = 0

            def getFileName(self):
                raise AssertionError(
                    "CISTEM streaming must use the logical Set API, "
                    "not reconstruct a Set from a persistence filename."
                )

            def loadAllProperties(self):
                self.reloads += 1

            def iterItems(self):
                return iter(self.items)

            def __iter__(self):
                return self.iterItems()

            def isStreamClosed(self):
                return self.closed

        class _LogicalPickingHarness(CistemProtFindParticles):
            def __init__(self, micSet, ctfSet):
                self.micDict = OrderedDict()
                self._micSet = micSet
                self.ctfRelations = _Pointer(ctfSet)

            def getInputMicrographs(self):
                return self._micSet

            def debug(self, *args, **kwargs):
                pass

        mic = _LogicalMic(1, "mic_001")
        micSet = _LogicalSet([mic])
        ctfSet = _LogicalSet([_LogicalCtf(mic)])
        protocol = _LogicalPickingHarness(micSet, ctfSet)

        readyMics, streamClosed = protocol._loadInputList()

        self.assertTrue(streamClosed)
        self.assertEqual(["mic_001"], list(readyMics))
        self.assertIsNotNone(readyMics["mic_001"].ctf)
        self.assertGreaterEqual(micSet.reloads, 1)
        self.assertGreaterEqual(ctfSet.reloads, 1)


if __name__ == "__main__":
    unittest.main()
