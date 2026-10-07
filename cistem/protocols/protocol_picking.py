# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
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
from collections import OrderedDict
from datetime import datetime

import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.constants import PROD
import pyworkflow.utils as pwutils
from pyworkflow.utils.properties import Message
import pwem.objects as emobj
from pwem.constants import RELATION_CTF
from pwem.protocols import ProtParticlePickingAuto
from pwem import emlib

from cistem import Plugin
from ..convert import readSetOfCoordinates, writeReferences
from ..constants import LOW_VARIANCE, FIND_PARTICLES_BIN
from .protocol_streaming_base import CistemStreamingBase

PICKING_STEP_NAMES = ('pickMicrographStep', 'pickMicrographListStep')


class CistemProtFindParticles(CistemStreamingBase, ProtParticlePickingAuto):
    """ Protocol to pick particles (ab-initio or reference-based) using cisTEM. """
    _label = 'find particles'
    _devStatus = PROD
    stepsExecutionMode = STEPS_PARALLEL

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        ProtParticlePickingAuto._defineParams(self, form)
        form.addParam('ctfRelations', params.RelationParam,
                      relationName=RELATION_CTF,
                      attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to the '
                           'input micrographs.')
        form.addParam('pickType', params.EnumParam, default=0,
                      important=True,
                      label='Picking algorithm',
                      choices=['Ab-initio', 'Reference-based'],
                      display=params.EnumParam.DISPLAY_HLIST)
        form.addParam('inputRefs', params.PointerParam,
                      condition='pickType==1', important=True,
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      label='Input references',
                      help='Provide a set of 2D templates to use in '
                           'the search.')
        form.addParam('maxradius', params.FloatParam, default=120.0,
                      label='Max particle radius (A)',
                      help='In Angstroms, the maximum radius of the '
                            'particles to be found. This also determines '
                            'the minimum distance between picks.')
        form.addParam('radius', params.FloatParam, default=80.0,
                      label='Characteristic particle radius (A)',
                      help='In Angstroms, the radius within which most '
                           'of the density is enclosed. The template '
                           'for picking is a soft-edge disc, where '
                           'the edge is 5 pixels wide and this '
                           'parameter defines the radius at which '
                           'the cosine-edge template reaches 0.5.')
        form.addParam('threshold', params.FloatParam, default=6.0,
                      label='Threshold peak height',
                      help='Particle coordinates will be defined as '
                           'the coordinates of any peak in the search '
                           'function which exceeds this threshold. '
                           'In numbers of standard deviations '
                           'above expected noise variations in '
                           'the scoring function. See Sigworth '
                           '(2004) for definition.')
        form.addParam('avoidHighVar', params.BooleanParam, default=False,
                      label='Avoid high variance areas',
                      help='Avoid areas with abnormally high local variance. '
                           'This can be effective in avoiding edges '
                           'of support films or contamination.')
        form.addParam('ptclWhite', params.BooleanParam, default=False,
                      label='Particles are white on a dark background?')

        form.addSection(label='Expert options')
        form.addParam('highRes', params.FloatParam, default=30.0,
                      label='Highest resolution used in picking (A)',
                      help='The template and micrograph will be resampled '
                           '(by Fourier cropping) to a pixel size of half '
                           'the resolution given here. Note that the '
                           'information in the corners of the Fourier '
                           'transforms (beyond the Nyquist frequency) '
                           'remains intact, so that there is some small '
                           'risk of bias beyond this resolution.')
        form.addParam('minDist', params.IntParam, default=0,
                      label='Minimum distance from edges (px)',
                      help='No particle shall be picked closer than '
                           'this distance from the edges of the micrograph. '
                           'In pixels.')
        form.addParam('useRadAvg', params.BooleanParam, default=True,
                      condition='pickType==1',
                      label='Use radial averages of templates',
                      help='Say yes if the templates should be '
                           'rotationally averaged')
        form.addParam('rotateRef', params.IntParam, default=0,
                      condition='pickType==1 and not useRadAvg',
                      label='Rotate each template this many times',
                      help='If > 0, each template image will be '
                           'rotated this number of times and the '
                           'micrograph will be searched for the rotated '
                           'template.')
        form.addParam('avoidLocMean', params.BooleanParam, default=True,
                      label='Avoid areas with abnormal local mean',
                      help='Avoid areas with abnormally low or high '
                           'local mean. This can be effective to avoid '
                           'picking from, e.g., contaminating ice crystals, '
                           'support film.')
        form.addParam('bgBoxes', params.IntParam, default=30,
                      label='Number of background boxes',
                      help='Number of background areas to use in estimating '
                           'the background spectrum. The larger the number '
                           'of boxes, the more accurate the estimate should '
                           'be, provided that none of the background boxes '
                           'contain any particles to be picked.')
        form.addParam('bgAlgo', params.EnumParam, default=LOW_VARIANCE,
                      choices=['Lowest variance', 'Variance near mode'],
                      display=params.EnumParam.DISPLAY_COMBO,
                      label='Algorithm to find background areas',
                      help='Testing so far suggests that areas of lowest '
                           'variance in experimental micrographs should be '
                           'used to estimate the background spectrum. '
                           'However, when using synthetic micrographs this '
                           'can lead to bias in the spectrum estimation '
                           'and the alternative (areas with local variances '
                           'near the mean of the distribution of local '
                           'variances) seems to perform better')

        self._defineStreamingParams(form)

        form.addParallelSection(threads=3)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self.inputStreaming = self.getInputMicrographs().isStreamOpen()

        if self.inputStreaming:
            # The conversion step must be scheduled here, not from inside
            # the generator: every picking step takes it as a prerequisite,
            # and a step inserted by a running step is not one the executor
            # has already planned around - the picking steps would wait on
            # it forever and the generator would poll with nothing to do.
            initialIds = self._insertInitialSteps()

            self._insertFunctionStep(
                self.resumableStepGeneratorStep,
                str(datetime.now()),
                prerequisites=initialIds,
                needsGPU=False,
            )
        else:
            # If not in streaming, then we will just insert a single step to
            # pick all micrographs at once since it is much faster
            self._insertInitialSteps()
            self._insertFunctionStep('_pickMicrographStep',
                                     self.getInputMicrographs(),
                                     *self._getPickArgs(),
                                     needsGPU=False)
            self._insertFunctionStep('createOutputStep', needsGPU=False)

            # Disable streaming functions:
            self._insertFinalSteps = self._doNothing
            self._stepsCheck = self._doNothing

    def resumableStepGeneratorStep(self, timestamp):
        """Run the streaming generator as a single resumable step."""
        self.stepsGeneratorStep()

    def stepsGeneratorStep(self):
        """Discover, pick and publish micrographs incrementally."""
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
        # The conversion step already ran before the generator started.
        self.initialIds = []

        self._restoreProcessedMicsFromPersistentState()

        while not self.finished:
            self._checkNewInput()
            self._checkNewOutput()

            if self.finished:
                break

            sleepOnWait = self._getStreamingSleepOnWait()
            if sleepOnWait > 0:
                self._streamingSleepOnWait()
            else:
                time.sleep(1)

    def _stepsCheck(self):
        """Persist steps created by the generator without legacy polling."""
        if getattr(self, '_newSteps', False):
            self.updateSteps()

    def _insertInitialSteps(self):
        """ Convert the input micrographs to mrc. """
        inputRefs = self.getInputReferences()
        refsId = inputRefs.strId() if inputRefs is not None else None
        convertId = self._insertFunctionStep('convertInputStep',
                                             self.getInputMicrographs().strId(),
                                             refsId, needsGPU=False)
        return [convertId]

    def _doNothing(self, *args):
        pass  # used to avoid some streaming functions

    def _loadSet(self, inputSet, SetClass, getKeyFunc, watermarkAttr,
                 knownIds):
        """Discover the items added to one input stream since the last poll.

        Only ids above that stream's watermark are queried and only those
        items are hydrated, so polling costs what just arrived rather than
        everything the stream has produced so far.
        """
        self.debug("Discovering new items from the logical input set.")

        newItems, producerClosed, terminalConsistent = (
            self._discoverNewInputItems(inputSet, watermarkAttr, knownIds))

        newItemDict = OrderedDict()

        for item in newItems:
            itemId = item.getObjId()

            if itemId in knownIds:
                continue

            knownIds.add(itemId)
            newItemDict[getKeyFunc(item)] = item

        return newItemDict, producerClosed and terminalConsistent

    def _loadMics(self, micSet):
        return self._loadSet(micSet, emobj.SetOfMicrographs,
                             lambda mic: mic.getMicName(),
                             '_lastMicId', self._knownMicIds)

    def _loadCTFs(self, ctfSet):
        return self._loadSet(ctfSet, emobj.SetOfCTF,
                             lambda ctf: ctf.getMicrograph().getMicName(),
                             '_lastCtfId', self._knownCtfIds)

    def _loadInputList(self):
        """Report the micrographs whose CTF has arrived.

        Both inputs are streams that advance independently, so each keeps
        its own watermark and whatever has no counterpart yet waits in a
        pending map. That map is the only thing re-examined per poll - a
        micrograph is never looked up in the input Set twice.
        """
        gapIds = getattr(self, '_resumeGapIds', None)

        if gapIds:
            micSet = self.getInputMicrographs()

            for mic in self._loadLogicalSetItemsByIds(micSet, gapIds):
                self._knownMicIds.add(mic.getObjId())
                self._micsWithoutCtf[mic.getMicName()] = mic

            self._resumeGapIds = set()

        newMics, micClosed = self._loadMics(self.getInputMicrographs())
        newCtfs, ctfClosed = self._loadCTFs(self.ctfRelations.get())

        self._micsWithoutCtf.update(newMics)

        for micKey, ctf in newCtfs.items():
            self._ctfByMicName[micKey] = ctf

        readyMics = OrderedDict()

        for micKey in list(self._micsWithoutCtf):
            ctf = self._ctfByMicName.pop(micKey, None)

            if ctf is None:
                continue

            mic = self._micsWithoutCtf.pop(micKey)
            mic.setCTF(ctf)
            readyMics[micKey] = mic

        if readyMics:
            scheduledNames = self._getScheduledPickingMicNames()

            for micKey in list(readyMics):
                if micKey in scheduledNames:
                    # Already has a step from an earlier run: it only needs
                    # publishing, so it must not be scheduled a second time.
                    self.micDict[micKey] = readyMics.pop(micKey)

        self._pendingMics.update(readyMics)

        return OrderedDict(self._pendingMics), micClosed and ctfClosed

    def _checkNewInput(self):
        """Discover ready micrographs directly from logical input Sets."""
        micDict, self.streamClosed = self._loadInputList()
        newMics = list(micDict.values())

        if newMics:
            self._insertNewMicsSteps(newMics)

            # pwem only takes whole batches; whatever it left out has to be
            # offered again, because discovery will not find it twice.
            for mic in newMics:
                if mic.getMicName() in self.micDict:
                    self._pendingMics.pop(mic.getMicName(), None)

            self.updateSteps()

    def _getFinishedPickingMicNames(self):
        """Return mic names represented by finished picking steps."""
        return self._collectStepArgKeys(PICKING_STEP_NAMES, keyType=str)

    def _getScheduledPickingMicNames(self):
        """Return mic names represented by persisted picking steps."""
        return self._collectStepArgKeys(PICKING_STEP_NAMES,
                                        onlyFinished=False, keyType=str)

    def _restoreProcessedMicsFromPersistentState(self):
        """Place the watermarks past what a previous run already handled.

        Coordinates carry the id of their micrograph, so the output Set
        says which micrographs were picked, and the step graph covers the
        ones that produced no coordinate at all. Discovery restarts above
        that, and anything below it that was never handled comes back as a
        gap rather than being walked for again on every poll.
        """
        self._lastMicId = getattr(self, '_lastMicId', 0)

        pickedMicIds = self._getPublishedPickingMicIds()

        if not pickedMicIds:
            return

        watermark, gapIds = self._resumeWatermarkWithGaps(
            self.getInputMicrographs(), pickedMicIds)

        self._lastMicId = max(self._lastMicId, watermark)
        self._knownMicIds.update(pickedMicIds)
        self._resumeGapIds = gapIds

    def _getPublishedPickingMicIds(self):
        """Micrograph ids already represented in the output coordinates.

        This is an id query on the output, not a walk over it, and it only
        runs once per execution.
        """
        micIds = self._getOutputUniqueValues(
            getattr(self, 'outputCoordinates', None), '_micId')

        if micIds is None:
            return set()

        return micIds

    def _checkNewOutput(self):
        """Publish finished picking steps without DONE sidecars."""
        if getattr(self, 'finished', False):
            return

        finishedNames = self._getFinishedPickingMicNames()

        # micDict holds what has been scheduled and not published yet, so
        # only that has to be looked at - never every micrograph seen.
        newDone = [
            mic for micName, mic in self.micDict.items()
            if micName in finishedNames
        ]

        allDone = (len(newDone) == len(self.micDict)
                   and not self._pendingMics
                   and not self._micsWithoutCtf)
        self.finished = self.streamClosed and allDone
        streamMode = pwobj.Set.STREAM_CLOSED if self.finished else pwobj.Set.STREAM_OPEN

        if newDone:
            self._updateOutputCoordSet(newDone, streamMode)

            for mic in newDone:
                self.micDict.pop(mic.getMicName(), None)
        elif not self.finished:
            if allDone:
                self._streamingSleepOnWait()
            return

        if self.finished:
            self._updateStreamState(streamMode)

    # --------------------------- STEPS functions -------------------------------
    def convertInputStep(self, micsId, refsId):
        """ Build the CTF lookup used as a fallback for non-streaming
        picking (those micrographs never go through _loadInputList, so
        they never get a CTF attached via mic.setCTF()), and convert
        whatever micrographs are already available.

        This step only runs once, near the start of the protocol, but
        a streaming input Set keeps growing afterwards - a micrograph
        (and its CTF) that arrives later would never be covered by this
        snapshot. _pickMicrographStep therefore also converts lazily,
        per micrograph, and prefers the always-fresh mic.getCTF()
        (kept current by _loadInputList on every poll) over this
        one-shot dict.
        """
        self.ctfDict = {}
        if self.ctfRelations.get() is not None:
            for ctf in self.ctfRelations.get():
                self.ctfDict[ctf.getMicrograph().getMicName()] = ctf.clone()

        for mic in self.getInputMicrographs():
            self._convertMic(mic)

        if refsId is not None:
            writeReferences(self.getInputReferences(),
                            self._getExtraPath('references.mrc'))

    def _convertMic(self, mic):
        """ Convert a micrograph to mrc in the tmp dir, if that has not
        already been done. Idempotent, so it is safe to call again for
        a micrograph convertInputStep already converted. """
        micName = mic.getFileName()
        outMic = os.path.join(self._getTmpPath(),
                              pwutils.replaceBaseExt(micName, 'mrc'))
        if not os.path.exists(outMic):
            if micName.endswith('.mrc'):
                pwutils.createAbsLink(os.path.abspath(micName), outMic)
            else:
                ih = emlib.image.ImageHandler()
                ih.convert(micName, outMic, emlib.DT_FLOAT)
        return outMic

    def _getMicCtf(self, mic):
        """ CTF for a micrograph reaching the picking step. Streaming
        input already has it attached (kept fresh by _loadInputList on
        every poll); non-streaming input never goes through
        _loadInputList, so fall back to the snapshot built once in
        convertInputStep. """
        ctf = mic.getCTF()
        if ctf is None:
            ctf = self.ctfDict.get(mic.getMicName())
        return ctf

    def pickMicrographStep(self, micName, *args):
        """Pick one micrograph without filesystem completion sidecars."""
        mic = self.micDict[micName]
        self.info("Picking micrograph: %s " % mic.getFileName())
        self._pickMicrograph(mic, *args)

    def pickMicrographListStep(self, micNameList, *args):
        """Pick a batch of micrographs without completion sidecars."""
        micList = [self.micDict[micName] for micName in micNameList]
        for mic in micList:
            self.info("Picking micrograph: %s " % mic.getFileName())
        self._pickMicrographList(micList, *args)

    def _pickMicrograph(self, mic, *args):
        self._pickMicrographStep([mic], *args)

    def _pickMicrographList(self, micList, *args):
        self._pickMicrographStep(micList, *args)

    def _pickMicrographStep(self, mics, args):
        """ Main func that runs the picking job.
        :param mics: micrograph list
        :param args: programs args
        """
        for mic in mics:
            try:
                outMic = self._convertMic(mic)
                ctf = self._getMicCtf(mic)
                if ctf is None:
                    raise Exception(
                        "No CTF available for micrograph %s"
                        % mic.getMicName())

                args.update({'micName': outMic,
                             'logFn': self._getLogFn(mic),
                             'outStack': self._getStackFn(mic),
                             'phaseShift': ctf.getPhaseShift() or 0.0,
                             'defocusU': ctf.getDefocusU(),
                             'defocusV': ctf.getDefocusV(),
                             'defocusAngle': ctf.getDefocusAngle()
                             })

                if self.pickType == 1:
                    args.update({
                        'refsFn': self._getExtraPath('references.mrc'),
                        'useRadAvg': 'YES' if self.useRadAvg else 'NO',
                        'rotateRef': self.rotateRef.get(),
                    })

                argsStr = self._getArgsStr()
                cmdArgs = argsStr % args

                self.runJob(self._getProgram(), cmdArgs,
                            env=Plugin.getEnviron())

                # Move output coords from tmp to extra
                pltFn = pwutils.replaceExt(self._getStackFn(mic), 'plt')
                pwutils.moveFile(pltFn, self._getPltFn(mic))

                # Clean tmp folder
                pwutils.cleanPath(outMic)
                pwutils.cleanPath(self._getLogFn(mic))
                pwutils.cleanPath(self._getStackFn(mic))
            except Exception as e:
                self.error("ERROR: Picking has failed for %s. %s" % (
                    mic.getMicName(), self._getErrorFromPickerTxt(mic, e)))
                self._writeFailedList([mic])

    def _getErrorFromPickerTxt(self, mic, e):
        """ Parse output log for errors.
        :param mic: input mic object
        :return: the error string
        """
        file = self._getLogFn(mic)
        try:
            with open(file, "r") as fh:
                for line in fh.readlines():
                    if line.startswith("Error"):
                        return line.replace("Error:", "")
        except OSError:
            pass
        return e

    def createOutputStep(self):
        """ Read the coordinates and define outputs."""
        micSet = self.getInputMicrographs()
        coordSet = self._createSetOfCoordinates(micSet)
        self.readCoordsFromMics(self._getExtraPath(), micSet,
                                coordSet)

        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(self.inputMicrographs, coordSet)

    # --------------------------- INFO functions --------------------------------
    def _validateStreamingThreads(self):
        inputMics = self.getInputMicrographs()

        if (inputMics is not None
                and inputMics.isStreamOpen()
                and self.numberOfThreads.get() < 3):
            return [
                'FindParticles streaming requires at least 3 threads.'
            ]

        return []

    def _validate(self):
        errors = ProtParticlePickingAuto._validate(self)
        errors.extend(self._validateStreamingThreads())
        return errors

    def _summary(self):
        summary = list()
        summary.append("Number of input micrographs: %d"
                       % self.getInputMicrographs().getSize())
        if self.getOutputsSize() > 0:
            summary.append("Number of particles picked: %d"
                           % self.getCoords().getSize())
            summary.append("Particle size: %d" % self.getCoords().getBoxSize())
            summary.append("Threshold: %0.2f" % self.threshold)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def _citations(self):
        return ['Sigworth2004']

    def _methods(self):
        methodsMsgs = []
        if self.getInputMicrographs() is not None:
            methodsMsgs.append("Input micrographs %s of size %d."
                               % (self.getObjectTag(self.getInputMicrographs()),
                                  self.getInputMicrographs().getSize()))

        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append('%s: User picked %d particles with a particle '
                               'size of %d and threshold %0.2f.'
                               % (self.getObjectTag(output), output.getSize(),
                                  output.getBoxSize(), self.threshold))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    # --------------------------- UTILS functions -------------------------------
    def _getProgram(self):
        """ Return program binary. """
        return Plugin.getProgram(FIND_PARTICLES_BIN)

    def _getPickArgs(self):
        """ Format arguments to call find_particles program. """
        inputMics = self.getInputMicrographs()
        acq = inputMics.getAcquisition()
        sampling = inputMics.getSamplingRate()

        self.argsDict = {'samplingRate': sampling,
                         'voltage': acq.getVoltage(),
                         'cs': acq.getSphericalAberration(),
                         'ampContrast': acq.getAmplitudeContrast(),
                         'templates': 'NO' if self.pickType == 0 else 'YES',
                         'radius': self.radius.get(),
                         'maxradius': self.maxradius.get(),
                         'highRes': self.highRes.get(),
                         'boxSize': int(self.maxradius.get() / sampling),
                         'minDist': self.minDist.get(),
                         'threshold': self.threshold.get(),
                         'avoidHighVar': 'YES' if self.avoidHighVar else 'NO',
                         'avoidLocMean': 'YES' if self.avoidLocMean else 'NO',
                         'algorithm': self.bgAlgo.get(),
                         'bgBoxes': self.bgBoxes.get(),
                         'ptclWhite': 'YES' if self.ptclWhite else 'NO'
                         }

        return [self.argsDict]

    def _getArgsStr(self):
        argsStr = """ << eof > %(logFn)s
%(micName)s
%(samplingRate)f
%(voltage)f
%(cs)f
%(ampContrast)f
%(phaseShift)f
%(defocusU)f
%(defocusV)f
%(defocusAngle)f"""
        if self.pickType == 0:
            argsStr += """
NO
%(radius)f
%(maxradius)f
%(highRes)f
%(outStack)s
%(boxSize)d
%(minDist)d
%(threshold)f
%(avoidHighVar)s
%(avoidLocMean)s
%(algorithm)d
%(bgBoxes)d
%(ptclWhite)s
eof"""
        else:  # ref-based picking
            argsStr += """
YES
%(refsFn)s
%(useRadAvg)s"""

            if self.useRadAvg:
                argsStr += """
%(maxradius)f
%(highRes)f
%(outStack)s
%(boxSize)d
%(minDist)d
%(threshold)f
%(avoidHighVar)s
%(avoidLocMean)s
%(algorithm)d
%(bgBoxes)d
%(ptclWhite)s
eof"""
            else:
                if self.rotateRef > 0:
                    argsStr += """
YES
%(rotateRef)d
%(maxradius)f
%(highRes)f
%(outStack)s
%(boxSize)d
%(minDist)d
%(threshold)f
%(avoidHighVar)s
%(avoidLocMean)s
%(algorithm)d
%(bgBoxes)d
%(ptclWhite)s
eof"""
                else:
                    argsStr += """
NO
%(maxradius)f
%(highRes)f
%(outStack)s
%(boxSize)d
%(minDist)d
%(threshold)f
%(avoidHighVar)s
%(avoidLocMean)s
%(algorithm)d
%(bgBoxes)d
%(ptclWhite)s
eof"""

        return argsStr

    def readCoordsFromMics(self, extraDir, micList, coordSet):
        """ Overwrite base class function. """
        if coordSet.getBoxSize() is None:
            coordSet.setBoxSize(128)

        readSetOfCoordinates(self._getExtraPath(), micList, coordSet)

    def _getMicrographDir(self, mic):
        """ Return an unique dir name for results of the micrograph. """
        return self._getTmpPath('mic_%06d' % mic.getObjId())

    def getInputMicrographs(self):
        return self.inputMicrographs.get()

    def _getLogFn(self, mic):
        """ Return output log file. """
        micName = mic.getFileName()
        return os.path.join(self._getTmpPath(),
                            pwutils.replaceBaseExt(micName, 'log'))

    def _getStackFn(self, mic):
        return self._getTmpPath('mic_%06d.mrc' % mic.getObjId())

    def _getPltFn(self, mic):
        """ Return output plt coords file. """
        micName = mic.getFileName()
        return os.path.join(self._getExtraPath(),
                            pwutils.replaceBaseExt(micName, 'plt'))

    def _getAllFailed(self):
        return self._getExtraPath('FAILED_all.TXT')

    def _writeFailedList(self, micList):
        """Do not persist failed micrographs in filesystem sidecars."""
        pass

    def getInputReferences(self):
        return self.inputRefs.get() if self.inputRefs.hasValue() else None
