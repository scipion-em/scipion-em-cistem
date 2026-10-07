# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es) [1]
# *              Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca) [2]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [3]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] Department of Anatomy and Cell Biology, McGill University
# * [3] MRC Laboratory of Molecular Biology, MRC-LMB
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
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import time
from datetime import datetime
from math import ceil
from threading import Thread

import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils

from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
from pyworkflow.gui.plotter import Plotter
from pwem.objects import Image
from pwem.protocols import ProtAlignMovies

from cistem import Plugin
from ..convert import readShiftsMovieAlignment
from ..constants import UNBLUR_BIN
from .protocol_streaming_base import CistemStreamingBase

# Module level: the output-id helper is called unbound on light test
# harnesses, so it cannot rely on a class attribute.
UNBLUR_OUTPUT_NAMES = (
    'outputMicrographs',
    'outputMicrographsDoseWeighted',
    'outputMovies',
    'outputMicrographsEven',
    'outputMicrographsOdd',
)


class CistemProtUnblur(CistemStreamingBase, ProtAlignMovies):
    """ This protocol wraps unblur movie alignment program. """

    _label = 'unblur'
    _devStatus = PROD
    CONVERT_TO_MRC = 'mrc'
    stepsExecutionMode = STEPS_PARALLEL

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif'] else 'mrc'

    def _defineAlignmentParams(self, form):
        form.addHidden('doSaveAveMic', params.BooleanParam,
                       default=True)
        form.addHidden('useAlignToSum', params.BooleanParam,
                       default=True)

        group = form.addGroup('Alignment')
        line = group.addLine('Frames to ALIGN',
                             help='Frames range to ALIGN on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to align, it means that you will '
                                  'align until the last frame of the movie.')
        line.addParam('alignFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('alignFrameN', params.IntParam, default=0,
                      label='to')

        group.addParam('binFactor', params.FloatParam, default=1.,
                       label='Binning factor',
                       help='1x or 2x. Bin stack before processing.')

        form.addParam('doComputePSD', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Compute PSD?",
                      help="If Yes, the protocol will compute for each "
                           "aligned micrograph the PSD using EMAN2.")
        form.addParam('doComputeMicThumbnail', params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=False,
                      label='Compute micrograph thumbnail?',
                      help='When using this option, we will compute a '
                           'micrograph thumbnail with EMAN2 and keep it with the '
                           'micrograph object for visualization purposes. ')
        form.addParam('extraProtocolParams', params.StringParam, default='',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Additional protocol parameters',
                      help="Here you can provide some extra parameters for the "
                           "protocol, not the underlying unblur program."
                           "You can provide many options separated by space. "
                           "\n\n*Options:* \n\n"
                           "--use_worker_thread \n"
                           " Use an extra thread to compute"
                           " PSD and thumbnail. This will allow requires "
                           "an extra CPU. ")

        form.addSection(label='Expert Options')
        line = form.addLine('Shifts (A): ',
                            help='Min and max shifts during alignment.\n\n'
                                 'The minimum shift can be applied '
                                 'during the initial refinement stage. '
                                 'Its purpose is to prevent images aligning '
                                 'to detector artifacts that may be '
                                 'reinforced in the initial sum which is '
                                 'used as the first reference. It is '
                                 'applied only during the first alignment '
                                 'round, and is ignored after that.\n'
                                 'The maximum shift can be applied in any '
                                 'single alignment round. Its purpose is '
                                 'to avoid alignment to spurious noise '
                                 'peaks by not considering unreasonably '
                                 'large shifts.  This limit is applied '
                                 'during every alignment round, but only'
                                 ' for that round, such that it can be '
                                 'exceeded over a number of successive rounds.')
        line.addParam('minShiftInitSearch', params.FloatParam, default='2.0',
                      label='Min shift')
        line.addParam('OutRadShiftLimit', params.FloatParam, default='40.0',
                      label='Max shift')

        group = form.addGroup('Exposure filter')
        group.addParam('doApplyDoseFilter', params.BooleanParam, default=True,
                       label='Exposure filter sums?',
                       help='If selected the resulting aligned movie sums '
                            'will be calculated using the exposure filter '
                            'as described in Grant and Grigorieff (2015). '
                            'Pre-exposure and dose per frame '
                            'should  be specified during movies import.')
        group.addParam('doRestoreNoisePwr', params.BooleanParam,
                       default=True,
                       label='Restore power? ',
                       help='If selected, and the exposure filter is used '
                            'to calculate the sum then the sum will be '
                            'high pass filtered to restore the noise '
                            'power. This is essentially the denominator '
                            'of Eq. 9 in Grant and Grigorieff (2015).')

        group = form.addGroup('Convergence')
        group.addParam('terminShiftThreshold', params.FloatParam,
                       default=1.0,
                       label='Termination threshold (A)',
                       help='The frames will be iteratively aligned '
                            'until either the maximum number of '
                            'iterations is reached, or if after an '
                            'alignment round every frame was shifted '
                            'by less than this threshold.')
        group.addParam('maximumNumberIterations', params.IntParam,
                       default=20,
                       label='Max iterations',
                       help='The maximum number of iterations that '
                            'can be run for the movie alignment. '
                            'If reached, the alignment will stop '
                            'and the current best values will be taken.')

        group = form.addGroup('Filter')
        group.addParam('bfactor', params.FloatParam,
                       default=1500.,
                       label='B-factor (A^2)',
                       help='This B-Factor is applied to the reference sum '
                            'prior to alignment. It is intended to low-pass '
                            'filter the images in order to prevent '
                            'alignment to spurious noise peaks and '
                            'detector artifacts.')

        line = group.addLine('Mask central cross?',
                             help='If selected, the Fourier transform of '
                                  'the reference will be masked by a cross '
                                  'centred on the origin of the transform. '
                                  'This is intended to reduce the influence '
                                  'of detector artifacts which often have '
                                  'considerable power along the central cross.')
        line.addParam('HWHoriFourMask', params.IntParam, default=1,
                      label='Horiz. mask (px)')
        line.addParam('HWVertFourMask', params.IntParam, default=1,
                      label='Vert. mask (px)')

        form.addParallelSection(threads=3, mpi=1)

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self.insertedDict = {}
        self.samplingRate = self.inputMovies.get().getSamplingRate()

        convertStepId = self._insertFunctionStep('_convertInputStep', prerequisites=[])
        self.convertCIStep = [convertStepId]

        self._insertFunctionStep(
            'resumableStepGeneratorStep',
            str(datetime.now()),
            prerequisites=self.convertCIStep,
            needsGPU=False
        )

    def resumableStepGeneratorStep(self, timestamp):
        self.stepsGeneratorStep()

    def stepsGeneratorStep(self):
        self.insertedDict = getattr(self, 'insertedDict', {})
        self.newDeps = []
        self.streamClosed = False
        self.finished = False

        self._restoreProcessedMoviesFromPersistentState()

        while not self.finished:
            # A failed step makes the executor stop and then join every
            # thread, this generator included: keep polling and the run
            # hangs for good with nothing left to do.
            # Returning, not breaking: the terminal steps below would
            # only add work the executor can never get to.
            if self._streamingMustStop():
                return

            self._checkNewInput()
            self._checkNewOutput()

            if self.finished:
                break

            sleepOnWait = self._getStreamingSleepOnWait()
            if sleepOnWait > 0:
                self._streamingSleepOnWait()
            else:
                time.sleep(1)

        # The final steps are inserted once the stream is known to be
        # exhausted, so nothing has to WAIT on an externally unlocked step.
        finalSteps = self._insertFinalSteps(self.newDeps)
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=finalSteps, needsGPU=False)

    def _stepsCheck(self):
        if getattr(self, '_newSteps', False):
            self.updateSteps()

    def createOutputStep(self):
        """Report failed movies against everything that was discovered.

        pwem's version measures the output against ``listOfMovies``, which
        here only holds what is still in flight - by the time this runs it
        is empty, so every failure would go unreported.
        """
        output = None

        for _, outputSet in self.iterOutputAttributes():
            output = outputSet
            break

        if output is None:
            return

        discovered = getattr(self, '_discoveredMovieCount', 0)

        if output.getSize() == 0 and discovered != 0:
            raise Exception("All movies failed, didn't create outputMicrographs."
                            "Please review movie processing steps above.")
        elif output.getSize() < discovered:
            self.warning(pwutils.yellowStr(
                "WARNING - Failed to align %d movies."
                % (discovered - output.getSize())))

    def _convertInputStep(self):
        """Convert correction images without creating DONE sidecars."""
        movies = self.inputMovies.get()
        movies.setGain(self._ProtProcessMovies__convertCorrectionImage(movies.getGain()))
        movies.setDark(self._ProtProcessMovies__convertCorrectionImage(movies.getDark()))

    def _loadInputList(self):
        """Discover the movies added since the last poll.

        Only the new ids are queried and only those movies are hydrated,
        so the cost of a poll follows what just arrived rather than
        everything the stream has produced so far.
        """
        movieSet = self.inputMovies.get()
        self.debug("Discovering new movies from the logical input set.")

        newMovies, producerClosed, terminalConsistent = (
            self._discoverNewInputItems(movieSet, '_lastInputId',
                                        self._knownMovieIds))

        gapIds = getattr(self, '_resumeGapIds', None)

        if gapIds:
            newMovies = (self._loadLogicalSetItemsByIds(movieSet, gapIds)
                         + newMovies)
            self._resumeGapIds = set()

        for movie in newMovies:
            movieId = movie.getObjId()

            if movieId in self._knownMovieIds:
                continue

            self._knownMovieIds.add(movieId)
            self._pendingMovies[movieId] = movie
            self._discoveredMovieCount += 1

        self.streamClosed = producerClosed and terminalConsistent

        # pwem's createOutputStep reports how many movies failed, which it
        # derives from listOfMovies; keep it holding what is still in
        # flight and track the discovered total separately.
        self.listOfMovies = list(self._pendingMovies.values())

    def _checkNewInput(self):
        """Discover new movies from the logical input Set."""
        self._loadInputList()

        newMovies = any(movie.getObjId() not in self.insertedDict
                        for movie in self.listOfMovies)

        if newMovies:
            dependencies = self._insertNewMoviesSteps(self.insertedDict, self.listOfMovies)
            self.newDeps.extend(dependencies)
            self.updateSteps()

    def _getPublishedMovieIds(self):
        """Return movie ids already present in logical outputs.

        This asks each output Set for its ids instead of walking it and
        hydrating every item, which matters because the outputs grow with
        the stream. It is only needed once, when restoring state.
        """
        publishedIds = set()

        for outputName in UNBLUR_OUTPUT_NAMES:
            outputSet = getattr(self, outputName, None)

            if outputSet is None:
                continue

            publishedIds.update(self._getOutputIdSet(outputSet))

        publishedIds.discard(None)

        return publishedIds

    def _getScheduledMovieIds(self):
        """Return movie ids represented by persisted processMovieStep steps."""
        return self._collectStepArgKeys(('processMovieStep',),
                                        onlyFinished=False,
                                        dictField='object.id')

    def _restoreProcessedMoviesFromPersistentState(self):
        """Restore already published or scheduled movies before discovery.

        This is the one place that reads the whole published state, and it
        runs once per execution rather than once per poll. The watermark
        starts past whatever is already published so a Continue does not
        walk the movies a previous run already dealt with.
        """
        self._knownMovieIds = getattr(self, '_knownMovieIds', set())
        self._pendingMovies = getattr(self, '_pendingMovies', {})
        self._lastInputId = getattr(self, '_lastInputId', 0)
        self._discoveredMovieCount = getattr(self, '_discoveredMovieCount', 0)

        publishedIds = self._getPublishedMovieIds()
        restoredIds = publishedIds | self._getScheduledMovieIds()

        self._publishedAnyOutput = bool(publishedIds)

        for movieId in restoredIds:
            self.insertedDict.setdefault(movieId, movieId)

        self._knownMovieIds.update(restoredIds)
        self._discoveredMovieCount = max(self._discoveredMovieCount,
                                         len(self._knownMovieIds))

        watermark, gapIds = self._resumeWatermarkWithGaps(
            self.inputMovies.get(), restoredIds)

        self._lastInputId = max(self._lastInputId, watermark)

        # A movie below the watermark that was never processed still has to
        # be picked up; discovery itself will not look that far back again.
        self._resumeGapIds = gapIds

    def _getFinishedMovieIds(self):
        """Return movie ids represented by finished processMovieStep steps."""
        return self._collectStepArgKeys(('processMovieStep',),
                                        dictField='object.id')

    def _checkNewOutput(self):
        """Publish finished movies without filesystem completion markers."""
        if getattr(self, 'finished', False):
            return

        # Only movies still in flight can become publishable, so the scan
        # follows what is pending instead of everything discovered so far.
        finishedIds = self._getFinishedMovieIds()
        newDone = [movie for movieId, movie in self._pendingMovies.items()
                   if movieId in finishedIds]

        self._firstTimeOutput = not self._publishedAnyOutput
        self.finished = (self.streamClosed
                         and len(newDone) == len(self._pendingMovies))
        streamMode = pwobj.Set.STREAM_CLOSED if self.finished else pwobj.Set.STREAM_OPEN

        if not newDone and not self.finished:
            return

        self._updateOutputSets(newDone, streamMode)

        for movie in newDone:
            self._pendingMovies.pop(movie.getObjId(), None)

        if newDone:
            self._publishedAnyOutput = True

        self.listOfMovies = list(self._pendingMovies.values())

    def processMovieStep(self, movieDict, hasAlignment):
        """Process one movie using the protocol step status as completion state."""
        import os

        import pwem.objects as emobj
        import pyworkflow.object as pwobj
        from pwem import emlib

        movie = emobj.Movie()
        movie.setAcquisition(emobj.Acquisition())

        if hasAlignment:
            movie.setAlignment(emobj.MovieAlignment())

        movie.setAttributesFromDict(movieDict, setBasic=True, ignoreMissing=True)

        movieFolder = self._getOutputMovieFolder(movie)
        movieFn = movie.getFileName()
        movieName = os.path.basename(movieFn)

        if self._filterMovie(movie):
            pwutils.makePath(movieFolder)
            pwutils.createAbsLink(os.path.abspath(movieFn),
                                  os.path.join(movieFolder, movieName))

            if movieName.endswith('bz2'):
                newMovieName = movieName.replace('.bz2', '')
                if not os.path.exists(newMovieName):
                    self.runJob('bzip2', '-d -f %s' % movieName, cwd=movieFolder)

            elif movieName.endswith('tbz'):
                newMovieName = movieName.replace('.tbz', '.mrc')
                if not os.path.exists(newMovieName):
                    self.runJob('tar', 'jxf %s' % movieName, cwd=movieFolder)

            elif movieName.endswith('.txt'):
                movieTxt = os.path.join(movieFolder, movieName)
                with open(movieTxt) as handle:
                    movieOrigin = os.path.basename(os.readlink(movieFn))
                    newMovieName = movieName.replace('.txt', '.mrcs')
                    imageHandler = emlib.image.ImageHandler()
                    for index, line in enumerate(handle):
                        if line.strip():
                            inputFrame = os.path.join(movieOrigin, line.strip())
                            imageHandler.convert(
                                inputFrame,
                                (index + 1, os.path.join(movieFolder, newMovieName))
                            )
            else:
                newMovieName = movieName

            convertExt = self._getConvertExtension(newMovieName)
            correctGain = self._doCorrectGain()

            if convertExt or correctGain:
                inputMovieFn = os.path.join(movieFolder, newMovieName)
                if inputMovieFn.endswith('.em'):
                    inputMovieFn += ':ems'

                if convertExt:
                    newMovieName = pwutils.replaceExt(newMovieName, convertExt)
                else:
                    newMovieName = '%s_corrected.%s' % os.path.splitext(newMovieName)

                outputMovieFn = os.path.join(movieFolder, newMovieName)

                if self._doCorrectGain():
                    self.info("Correcting gain and dark '%s' -> '%s'"
                              % (inputMovieFn, outputMovieFn))
                    gain, dark = self.getGainAndDark()
                    self.correctGain(inputMovieFn, outputMovieFn,
                                     gainFn=gain, darkFn=dark)
                else:
                    self.info("Converting movie '%s' -> '%s'"
                              % (inputMovieFn, outputMovieFn))
                    emlib.image.ImageHandler().convertStack(inputMovieFn, outputMovieFn)

            movie._originalFileName = pwobj.String(objDoStore=False)
            movie._originalFileName.set(movie.getFileName())
            movie.setFileName(os.path.join(movieFolder, newMovieName))
            self.info("Processing movie: %s" % movie.getFileName())

            self._processMovie(movie)

            if self._doMovieFolderCleanUp():
                self._cleanMovieFolder(movieFolder)

    def _cleanMovieFolder(self, movieFolder):
        """Remove a movie's working folder without going through a shell.

        pwem builds a shell command out of the folder path, so a project
        path with a space in it turns the cleanup into the deletion of
        whatever the shell reads as a second argument. Refuse anything
        that is not inside this run's own working directory, and remove
        it through the filesystem API rather than a command line.
        """
        if pwutils.envVarOn('SCIPION_DEBUG_NOCLEAN'):
            self.info('Clean movie data DISABLED. '
                      'Movie folder will remain in disk!!!')
            return

        workspace = os.path.realpath(self._getTmpPath())
        target = os.path.realpath(movieFolder)

        if target != workspace and not target.startswith(workspace + os.sep):
            self.warning("Refusing to remove %s: it is outside this run's "
                         "working directory." % movieFolder)
            return

        if target == workspace:
            self.warning("Refusing to remove the working directory itself: "
                         "every other movie in flight lives there too.")
            return

        self.info("Erasing movie folder: %s" % movieFolder)
        pwutils.cleanPath(target)

    def _processMovie(self, movie):
        inputMovies = self.getInputMovies()

        try:
            # processMovieStep (the pwem base class step calling this
            # hook) has no exception boundary of its own around it - a
            # single corrupted/unusual movie failing here (e.g. a bad
            # tiff link or an acquisition attribute missing while
            # building the unblur arguments) must not crash the whole
            # protocol, so these are inside the try like every other
            # failure in this function.
            self._createTifLink(movie)
            self._argsUnblur(movie)

            self.runJob(self._getProgram(), self._args, env=Plugin.getEnviron())

            def _extraWork():
                outMicFn = self._getMicFn(movie)

                if self.doComputePSD:
                    self._computePSD(outMicFn, outputFn=self._getPsdCorr(movie))

                self._saveAlignmentPlots(movie, inputMovies.getSamplingRate())

                if self._doComputeMicThumbnail():
                    self.computeThumbnail(outMicFn,
                                          outputFn=self._getOutputMicThumbnail(movie))

            if self._useWorkerThread():
                # The side thread is joined before this step returns. What
                # a FINISHED step says is that the movie is done, and
                # _checkNewOutput publishes it on that word alone - a
                # thread still writing the PSD, the plots or the thumbnail
                # would have the micrograph published pointing at files
                # that are not there yet. Movies already overlap with each
                # other: each one is a step of its own.
                thread = Thread(target=_extraWork)
                thread.start()
                self._trackWorkerThread(thread)
                thread.join()
            else:
                _extraWork()

        except Exception as e:
            self.error("ERROR: Unblur has failed for %s. %s" % (
                self._getMovieFn(movie), self._getErrorFromUnblurTxt(movie, e)))

    def _getErrorFromUnblurTxt(self, movie, e):
        """ Parse output log for errors.
        :param movie: input movie object
        :return: the error string
        """
        file = self._getShiftsFn(movie)
        try:
            with open(file, "r") as fh:
                for line in fh.readlines():
                    if line.startswith("Error"):
                        return line.replace("Error:", "")
        except OSError:
            pass
        return e

    def _insertFinalSteps(self, deps):
        stepId = self._insertFunctionStep('waitForThreadStep',
                                          prerequisites=deps, needsGPU=False)
        return [stepId]

    def _trackWorkerThread(self, thread):
        """Remember a side thread so nothing can outlive the run.

        Movie steps run in parallel, so the registry is guarded by the
        lock every Protocol already owns - a lock of its own would have to
        live on the protocol, which is not somewhere a thread primitive
        belongs.
        """
        with self._lock:
            live = [
                tracked for tracked in getattr(self, '_workerThreads', [])
                if tracked.is_alive()
            ]
            live.append(thread)
            self._workerThreads = live

    def waitForThreadStep(self):
        """Join whatever side threads are still running.

        Each movie's step already joins its own thread, so by the time
        this runs there is normally nothing left. It stays as the explicit
        guarantee that the run does not end with work still in flight -
        the previous version slept for a fixed ten seconds instead, which
        was neither a guarantee nor free.
        """
        with self._lock:
            threads = list(getattr(self, '_workerThreads', []))

        for thread in threads:
            thread.join()

    # --------------------------- INFO functions -------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, 'outputMicrographs') or \
                hasattr(self, 'outputMicrographsDoseWeighted'):
            summary.append('Aligned %d movies using unblur.'
                           % self.inputMovies.get().getSize())
        else:
            summary.append('Output is not ready')

        return summary

    def _citations(self):
        return ["Campbell2012", "Grant2015b"]

    def _validateStreamingThreads(self):
        if self.numberOfThreads.get() < 3:
            return ['Unblur streaming requires at least 3 threads.']
        return []

    def _validate(self):
        # Check base validation before the specific ones
        errors = ProtAlignMovies._validate(self)
        errors.extend(self._validateStreamingThreads())

        if self.doApplyDoseFilter and self.inputMovies.get():
            inputMovies = self.inputMovies.get()
            doseFrame = inputMovies.getAcquisition().getDosePerFrame()

            if doseFrame == 0.0 or doseFrame is None:
                errors.append('Dose per frame for input movies is 0 or not '
                              'set. You cannot apply dose filter.')

        if self.doComputeMicThumbnail or self.doComputePSD:
            try:
                from pwem import Domain
                _ = Domain.importFromPlugin('eman2', doRaise=True)
            except:
                errors.append("EMAN2 plugin not found!\nComputing thumbnails "
                              "or PSD requires EMAN2 plugin and binaries installed.")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getProgram(self):
        """ Return program binary. """
        return Plugin.getProgram(UNBLUR_BIN)

    def _argsUnblur(self, movie):
        """ Format arguments to call unblur program. """
        inputMovies = self.getInputMovies()
        if self.doApplyDoseFilter:
            preExp, dose = self._getCorrectedDose(inputMovies)
        else:
            preExp, dose = 0.0, 0.0

        args = {'movieName': self._getMovieFn(movie),
                'micFnName': self._getMicFn(movie),
                'shiftsFn': self._getShiftsFn(movie),
                'samplingRate': self.samplingRate,
                'voltage': movie.getAcquisition().getVoltage(),
                'bfactor': self.bfactor.get(),
                'minShiftInitSearch': self.minShiftInitSearch.get(),
                'OutRadShiftLimit': self.OutRadShiftLimit.get(),
                'HWVertFourMask': self.HWVertFourMask.get(),
                'HWHoriFourMask': self.HWHoriFourMask.get(),
                'terminShiftThreshold': self.terminShiftThreshold.get(),
                'maximumNumberIterations': self.maximumNumberIterations.get(),
                'applyDoseFilter': 'YES' if self.doApplyDoseFilter else 'NO',
                'doRestoreNoisePwr': 'YES' if self.doRestoreNoisePwr else 'NO',
                'exposurePerFrame': dose,
                'binFactor': self.binFactor.get(),
                'alignFrame0': self.alignFrame0.get(),
                'alignFrameN': self.alignFrameN.get(),
                'gainCorrected': 'NO' if inputMovies.getGain() else 'YES',
                'gainFn': inputMovies.getGain(),
                'preExposureAmount': preExp
                }

        argsStr = """ << eof > %(shiftsFn)s
%(movieName)s
%(micFnName)s
%(samplingRate)f
%(binFactor)f
%(applyDoseFilter)s"""

        if self.doApplyDoseFilter:
            argsStr += """
%(voltage)f
%(exposurePerFrame)f
%(preExposureAmount)f"""

        argsStr += """
YES
%(minShiftInitSearch)f
%(OutRadShiftLimit)f
%(bfactor)f
%(HWVertFourMask)d
%(HWHoriFourMask)d
%(terminShiftThreshold)f
%(maximumNumberIterations)d"""

        if self.doApplyDoseFilter:
            argsStr += """
%(doRestoreNoisePwr)s"""

        if inputMovies.getGain():
            argsStr += """
%(gainCorrected)s
%(gainFn)s
%(alignFrame0)d
%(alignFrameN)d
NO
eof\n
"""
        else:
            argsStr += """
%(gainCorrected)s
%(alignFrame0)d
%(alignFrameN)d
NO
eof\n
"""

        self._args = argsStr % args

    def _getMovieFn(self, movie):
        movieFn = movie.getFileName()
        if movieFn.endswith("tiff"):
            return pwutils.replaceExt(movieFn, "tif")
        else:
            return movieFn

    def _createTifLink(self, movie):
        # unblur recognises only tif, not tiff
        movieFn = movie.getFileName()
        if movieFn.endswith("tiff"):
            pwutils.createLink(movieFn, self._getMovieFn(movie))

    def _getMicFn(self, movie):
        if self.doApplyDoseFilter:
            return self._getExtraPath(self._getOutputMicWtName(movie))
        else:
            return self._getExtraPath(self._getOutputMicName(movie))

    def _getShiftsFn(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_shifts.txt')

    def _getMovieShifts(self, movie):
        """ Returns the x and y shifts for the alignment of this movie. """
        pixSize = movie.getSamplingRate()
        shiftFn = self._getShiftsFn(movie)
        xShifts, yShifts = readShiftsMovieAlignment(shiftFn)
        # convert shifts from Angstroms to px
        xShiftsCorr = [x / pixSize for x in xShifts]
        yShiftsCorr = [y / pixSize for y in yShifts]

        return xShiftsCorr, yShiftsCorr

    def _doComputeMicThumbnail(self):
        return self.doComputeMicThumbnail

    def _computePSD(self, inputFn, outputFn, scaleFactor=6):
        """ Generate a thumbnail of the PSD with EMAN2"""
        args = "%s %s " % (inputFn, outputFn)
        args += "--process=math.realtofft --meanshrink %s " % scaleFactor
        args += "--fixintscaling=sane"

        from pwem import Domain
        eman2 = Domain.importFromPlugin('eman2')
        from pyworkflow.utils.process import runJob
        runJob(self._log, eman2.Plugin.getProgram('e2proc2d.py'), args,
               env=eman2.Plugin.getEnviron())

        return outputFn

    def _preprocessOutputMicrograph(self, mic, movie):
        mic.plotGlobal = Image(location=self._getPlotGlobal(movie))
        if self.doComputePSD:
            mic.psdCorr = Image(location=self._getPsdCorr(movie))
        if self._doComputeMicThumbnail():
            mic.thumbnail = Image(location=self._getOutputMicThumbnail(movie))

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _getPlotGlobal(self, movie):
        return self._getNameExt(movie, '_global_shifts', 'png', extra=True)

    def _getPsdCorr(self, movie):
        return self._getNameExt(movie, '_psd', 'png', extra=True)

    def _saveAlignmentPlots(self, movie, pixSize):
        """ Compute alignment shift plots and save to file as png images. """
        shiftsX, shiftsY = self._getMovieShifts(movie)
        first, _ = self._getFrameRange(movie.getNumberOfFrames(), 'align')
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first, pixSize)
        plotter.savefig(self._getPlotGlobal(movie))
        plotter.close()

    def _useWorkerThread(self):
        return '--use_worker_thread' in self.extraProtocolParams.get()

    def getInputMovies(self):
        return self.inputMovies.get()

    def _createOutputMicrographs(self):
        return not self.doApplyDoseFilter

    def _createOutputWeightedMicrographs(self):
        return self.doApplyDoseFilter


def createGlobalAlignmentPlot(meanX, meanY, first, pixSize):
    """ Create a plotter with the shift per frame. """
    sumMeanX = []
    sumMeanY = []

    def px_to_ang(apx):
        y1, y2 = apx.get_ylim()
        x1, x2 = apx.get_xlim()
        ax_ang2.set_ylim(y1*pixSize, y2*pixSize)
        ax_ang.set_xlim(x1*pixSize, x2*pixSize)
        ax_ang.figure.canvas.draw()
        ax_ang2.figure.canvas.draw()

    figureSize = (6, 4)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()
    ax_px = figure.add_subplot(111)
    ax_px.grid()
    ax_px.set_xlabel('Shift x (px)')
    ax_px.set_ylabel('Shift y (px)')

    ax_ang = ax_px.twiny()
    ax_ang.set_xlabel('Shift x (A)')
    ax_ang2 = ax_px.twinx()
    ax_ang2.set_ylabel('Shift y (A)')

    i = first
    skipLabels = ceil(len(meanX)/10.0)
    labelTick = 1

    for x, y in zip(meanX, meanY):
        sumMeanX.append(x)
        sumMeanY.append(y)
        if labelTick == 1:
            ax_px.text(x - 0.02, y + 0.02, str(i))
            labelTick = skipLabels
        else:
            labelTick -= 1
        i += 1

    # automatically update lim of ax_ang when lim of ax_px changes.
    ax_px.callbacks.connect("ylim_changed", px_to_ang)
    ax_px.callbacks.connect("xlim_changed", px_to_ang)

    ax_px.plot(sumMeanX, sumMeanY, color='b')
    ax_px.plot(sumMeanX, sumMeanY, 'yo')
    ax_px.plot(sumMeanX[0], sumMeanY[0], 'ro', markersize=10, linewidth=0.5)
    ax_px.set_title('Global frame alignment')

    plotter.tightLayout()

    return plotter
