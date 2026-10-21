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
"""A movie's step must not finish before the movie's artefacts exist.

Unblur can compute the PSD, the alignment plots and the thumbnail on a
side thread. The step that owns the movie returns without waiting for
it, so the step is FINISHED - and _checkNewOutput publishes the
micrograph - while those files are still being written. Whatever the
step reports as done has to be done.
"""
import threading
import time
import unittest

from cistem.protocols.protocol_picking import CistemProtFindParticles
from cistem.protocols.protocol_unblur import CistemProtUnblur


class _Movie:
    def __init__(self, objId=1):
        self._objId = objId

    def getObjId(self):
        return self._objId

    def getFileName(self):
        return 'movie_%06d.mrc' % self._objId


class _InputMovies:
    def getSamplingRate(self):
        return 1.0


class _ProcessHarness:
    """Drives the real _processMovie with the slow work stubbed out."""

    _processMovie = CistemProtUnblur._processMovie
    waitForThreadStep = CistemProtUnblur.waitForThreadStep
    _insertFinalSteps = CistemProtUnblur._insertFinalSteps
    _trackWorkerThread = CistemProtUnblur._trackWorkerThread

    def __init__(self, useWorkerThread=True, extraWorkSeconds=0.2):
        # The real protocol guards the registry with the lock every
        # Protocol already owns.
        self._lock = threading.RLock()
        self._workerThreads = []
        self._useWorker = useWorkerThread
        self._extraWorkSeconds = extraWorkSeconds
        self.doComputePSD = True
        self._args = ''
        self.finishedArtefacts = []
        self.errors = []
        self.insertedSteps = []

    # -- the pieces _processMovie leans on -------------------------------
    def getInputMovies(self):
        return _InputMovies()

    def _createTifLink(self, movie):
        pass

    def _argsUnblur(self, movie):
        pass

    def _getProgram(self):
        return 'unblur'

    def runJob(self, *args, **kwargs):
        pass

    def _getMicFn(self, movie):
        return 'mic_%06d.mrc' % movie.getObjId()

    def _getPsdCorr(self, movie):
        return 'psd_%06d.mrc' % movie.getObjId()

    def _computePSD(self, inputFn, outputFn):
        time.sleep(self._extraWorkSeconds)
        self.finishedArtefacts.append(outputFn)

    def _saveAlignmentPlots(self, movie, samplingRate):
        self.finishedArtefacts.append('plots_%06d' % movie.getObjId())

    def _doComputeMicThumbnail(self):
        return True

    def computeThumbnail(self, inputFn, outputFn):
        self.finishedArtefacts.append(outputFn)

    def _getOutputMicThumbnail(self, movie):
        return 'thumb_%06d.png' % movie.getObjId()

    def _useWorkerThread(self):
        return self._useWorker

    def _getMovieFn(self, movie):
        return movie.getFileName()

    def _getErrorFromUnblurTxt(self, movie, error):
        return str(error)

    def error(self, message):
        self.errors.append(message)

    def _insertFunctionStep(self, *args, **kwargs):
        self.insertedSteps.append(args[0])
        return len(self.insertedSteps)


class TestTheStepOwnsItsArtefacts(unittest.TestCase):

    def testArtefactsAreCompleteWhenTheStepReturns(self):
        harness = _ProcessHarness(useWorkerThread=True)

        harness._processMovie(_Movie(1))

        self.assertEqual(
            len(harness.finishedArtefacts),
            3,
            "The step returned while its side thread was still writing "
            "the PSD, the plots and the thumbnail: the movie gets "
            "published pointing at files that are not there yet.",
        )
        self.assertEqual(harness.errors, [])

    def testNoSideThreadOutlivesTheStep(self):
        before = threading.active_count()
        harness = _ProcessHarness(useWorkerThread=True)

        harness._processMovie(_Movie(1))

        self.assertEqual(
            threading.active_count(),
            before,
            "A side thread left running past the end of the step keeps "
            "writing into a movie the protocol already considers done.",
        )

    def testTheSynchronousPathIsUnchanged(self):
        harness = _ProcessHarness(useWorkerThread=False)

        harness._processMovie(_Movie(1))

        self.assertEqual(len(harness.finishedArtefacts), 3)

    def testAFailingMovieIsStillIsolated(self):
        harness = _ProcessHarness(useWorkerThread=True)
        harness._argsUnblur = self._raise

        harness._processMovie(_Movie(1))

        self.assertEqual(
            len(harness.errors),
            1,
            "One bad movie must be reported, not raised out of the step.",
        )

    @staticmethod
    def _raise(movie):
        raise RuntimeError('argument preparation failed')


class TestTheFinalWaitIsNotABlindSleep(unittest.TestCase):

    def testWaitingForThreadsDoesNotCostAFixedDelay(self):
        harness = _ProcessHarness(useWorkerThread=True)

        started = time.time()
        harness.waitForThreadStep()
        elapsed = time.time() - started

        self.assertLess(
            elapsed,
            5,
            "The terminal wait sleeps for a fixed ten seconds whether or "
            "not anything is outstanding, which is neither a guarantee "
            "nor free.",
        )


if __name__ == '__main__':
    unittest.main()


class _ConvertHarness:
    """Drives FindParticles' one-shot conversion step."""

    convertInputStep = CistemProtFindParticles.convertInputStep

    def __init__(self, inputStreaming):
        self.inputStreaming = inputStreaming
        self.converted = []
        self.ctfWalks = 0

    def getInputMicrographs(self):
        return []

    def _convertMic(self, mic):
        self.converted.append(mic)

    def getInputReferences(self):
        return None

    class _CtfRelations:
        def __init__(self, owner):
            self._owner = owner

        def get(self):
            return self

        def __iter__(self):
            self._owner.ctfWalks += 1
            return iter(())

    @property
    def ctfRelations(self):
        return _ConvertHarness._CtfRelations(self)


class TestConversionDoesNotBuildDeadState(unittest.TestCase):
    """The CTF lookup is only ever read on the non-streaming path.

    While streaming, every micrograph reaching a picking step already
    carries the CTF that _loadInputList attached to it, so building the
    lookup walks and clones the whole CTF Set for something nothing
    reads.
    """

    def testStreamingDoesNotWalkTheCtfRelations(self):
        harness = _ConvertHarness(inputStreaming=True)

        harness.convertInputStep('mics', None)

        self.assertEqual(
            harness.ctfWalks,
            0,
            "While streaming the lookup is never read: building it only "
            "hydrates and retains a clone of every CTF for nothing.",
        )

    def testTheNonStreamingPathStillBuildsIt(self):
        harness = _ConvertHarness(inputStreaming=False)

        harness.convertInputStep('mics', None)

        self.assertEqual(
            harness.ctfWalks,
            1,
            "Micrographs picked outside streaming never go through "
            "_loadInputList, so the lookup is the only CTF they get.",
        )
