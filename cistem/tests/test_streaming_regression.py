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
    def __init__(self, obj_id=1):
        self._obj_id = obj_id

    def getObjId(self):
        return self._obj_id

    def getFileName(self):
        return "mic_001.mrc"

    def getMicName(self):
        return "mic_001"


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
