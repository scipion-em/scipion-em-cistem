# **************************************************************************
# *
# * Authors:     Josue Gomez BLanco (josue.gomez-blanco@mcgill.ca) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [3]
# *
# * [1] Department of Anatomy and Cell Biology, McGill University
# * [2] SciLifeLab, Stockholm University
# * [3] MRC Laboratory of Molecular Biology (MRC-LMB)
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

import json
import os
import time
from collections import OrderedDict
from datetime import datetime

import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD
from pyworkflow.object import Boolean, Set
from pwem.protocols import ProtCTFMicrographs
from pwem.objects import CTFModel
from pwem import emlib

from .program_ctffind import ProgramCtffind


class CistemProtCTFFind(ProtCTFMicrographs):
    """ Estimate CTF for a set of micrographs with ctffind.
    
    To find more information about ctffind visit:
    https://grigoriefflab.umassmed.edu/ctffind4
    """
    _label = 'ctffind'
    _devStatus = PROD
    recalculate = Boolean(False, objDoStore=False)  # Legacy July 2024: to fake old recalculate param
    # that is still used in the ProtCTFMicrographs (to be removed)

    def _defineParams(self, form):
        ProgramCtffind.defineInputParams(form)
        ProgramCtffind.defineProcessParams(form)
        self._defineStreamingParams(form)

    def _defineCtfParamsDict(self):
        ProtCTFMicrographs._defineCtfParamsDict(self)
        self._ctfProgram = ProgramCtffind(self)

    # -------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        """Insert only the resumable streaming generator."""
        self._insertFunctionStep(self.resumableStepGeneratorStep,
                                 str(datetime.now()),
                                 needsGPU=False)

    def resumableStepGeneratorStep(
            self,
            timestamp,
    ):
        """Run the generator as a unique step on every resume."""
        self.stepsGeneratorStep()

    def stepsGeneratorStep(self):
        """Discover, process and publish CTFs incrementally."""
        self._defineCtfParamsDict()

        self.micDict = OrderedDict()
        self.streamClosed = False
        self.finished = False
        self.initialIds = self._insertInitialSteps()

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
                # The legacy streaming default is zero because polling was
                # driven externally by _stepsCheck. A generator must yield
                # CPU while waiting for new input or processing completion.
                time.sleep(1)

    def _stepsCheck(self):
        """Persist steps created by the generator, without legacy polling."""
        if getattr(
                self,
                "_newSteps",
                False,
        ):
            self.updateSteps()

    def _insertInitialSteps(self):
        """No filesystem initialization is required for streaming state."""
        return []

    def estimateCtfStep(self, micName, *args):
        """Estimate one CTF without filesystem completion markers.

        The persisted protocol-step status is the completion authority.
        """
        mic = self.micDict[micName]

        self.info("Estimating CTF of micrograph: %s " % mic.getObjId())

        self._estimateCTF(mic, *args)

    def estimateCtfListStep(
            self,
            micNameList,
            *args,
    ):
        """Estimate a CTF batch without filesystem completion markers."""
        micList = [self.micDict[micName] for micName in micNameList]

        self.info("Estimating CTF for micrographs: %s"
                  % [mic.getObjId() for mic in micList])

        self._estimateCtfList(micList, *args)

    def _doCtfEstimation(self, mic, **kwargs):
        """ Run ctffind with required parameters.
        :param mic: input mic object
        :param kwargs: dict with arguments
        """
        if self.usePowerSpectra:
            micFn = mic._powerSpectra.getFileName()
            powerSpectraPix = self.psSampling
        else:
            micFn = mic.getFileName()
            powerSpectraPix = None
        micDir = self._getTmpPath('mic_%06d' % mic.getObjId())
        micFnMrc = os.path.join(micDir, pwutils.replaceBaseExt(micFn, 'mrc'))

        try:
            # Create micrograph dir and convert here (instead of
            # outside this try) so a missing/corrupted micrograph is
            # caught and logged per-micrograph, like every other
            # failure in this function, instead of raising uncaught
            # and crashing the whole protocol via the single-mic
            # streaming step, which has no exception boundary of its
            # own around _estimateCTF.
            pwutils.makePath(micDir)

            if not os.path.exists(micFn):
                raise FileNotFoundError(
                    "Missing input micrograph: %s" % micFn)

            if micFn.endswith('.mrc'):
                pwutils.createAbsLink(os.path.abspath(micFn), micFnMrc)
            else:
                ih = emlib.image.ImageHandler()
                ih.convert(micFn, micFnMrc, emlib.DT_FLOAT)

            program, args = self._ctfProgram.getCommand(
                micFn=micFnMrc,
                powerSpectraPix=powerSpectraPix,
                ctffindOut=self._getCtfOutPath(mic),
                ctffindPSD=self._getPsdPath(mic),
                **kwargs)
            self.runJob(program, args)

            pwutils.cleanPath(micDir)

        except Exception as e:
            self.error("ERROR: Ctffind has failed for %s. %s" % (
                micFnMrc, self._getErrorFromCtffindTxt(mic, e)))

    def _getErrorFromCtffindTxt(self, mic, e):
        """ Parse output log for errors.
        :param mic: input mic object
        :return: the error string
        """
        file = self._getCtfOutPath(mic)
        try:
            with open(file, "r") as fh:
                for line in fh.readlines():
                    if "Error:" in line:
                        return line.split("Error:")[-1]
        except OSError:
            pass
        return e

    def _estimateCTF(self, mic, *args):
        """ Redefined func from the base class. """
        self._doCtfEstimation(mic)

    def _getPublishedCtfMicNames(self):
        """Return micrographs already present in the logical CTF output."""
        outputCtf = getattr(
            self,
            "outputCTF",
            None,
        )

        if outputCtf is None:
            return set()

        iterator = outputCtf.iterItems() if hasattr(outputCtf, "iterItems") else iter(outputCtf)

        publishedMicNames = set()

        for ctf in iterator:
            mic = ctf.getMicrograph()

            if mic is not None:
                publishedMicNames.add(mic.getMicName())

        return publishedMicNames

    def _getFinishedCtfMicNames(self):
        """Return micrographs represented by finished CTF processing steps."""
        finishedMicNames = set()

        for step in getattr(
                self,
                "_steps",
                [],
        ):
            if not step.isFinished():
                continue

            funcName = getattr(step, "funcName", None)

            if hasattr(funcName, "get"):
                funcName = funcName.get()

            if funcName not in (
                    "estimateCtfStep",
                    "estimateCtfListStep",
            ):
                continue

            argsStr = getattr(step, "argsStr", None)

            if hasattr(argsStr, "get"):
                argsStr = argsStr.get("[]")

            try:
                args = json.loads(argsStr or "[]")
            except (
                    TypeError,
                    ValueError,
            ):
                continue

            if not args:
                continue

            micArg = args[0]

            if (
                    funcName
                    == "estimateCtfListStep"
            ):
                if isinstance(micArg, list):
                    finishedMicNames.update(micArg)
            else:
                finishedMicNames.add(micArg)

        return finishedMicNames

    def _restoreProcessedMicsFromPersistentState(self):
        """Restore already processed inputs when resuming the generator.

        Published CTFs and finished processing steps are persistent state.
        Seeding micDict with their corresponding logical input objects keeps
        discovery incremental after Continue without filesystem markers.
        """
        processedMicNames = (self._getPublishedCtfMicNames()
                             | self._getFinishedCtfMicNames())

        if not processedMicNames:
            return

        inputMics = self.getInputMicrographs()
        inputMics.loadAllProperties()

        for mic in inputMics.iterItems():
            micName = mic.getMicName()

            if micName in processedMicNames:
                self.micDict[micName] = mic.clone()

    def _checkNewOutput(self):
        """Publish finished CTF steps without filesystem sidecars.

        Processing completion comes from persisted step state. Publication
        state comes from the logical output Set. This makes resume and
        streaming independent of DONE marker files.
        """
        if getattr(
            self,
            "finished",
            False,
        ):
            return

        publishedMicNames = self._getPublishedCtfMicNames()

        finishedMicNames = self._getFinishedCtfMicNames()

        listOfMics = list(self.micDict.values())

        newDone = [
            mic
            for mic in listOfMics
            if (
                mic.getMicName()
                in finishedMicNames
                and mic.getMicName()
                not in publishedMicNames
            )
        ]

        completedMicNames = (
            publishedMicNames
            | finishedMicNames
        )

        allDone = all(
            mic.getMicName()
            in completedMicNames
            for mic in listOfMics
        )

        self.finished = self.streamClosed and allDone

        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            self._updateOutputCTFSet(newDone, streamMode)
        elif not self.finished:
            if allDone:
                self._streamingSleepOnWait()

            return

        if self.finished:
            self._updateStreamState(streamMode)

            outputStep = self._getFirstJoinStep()

            if outputStep and outputStep.isWaiting():
                from pyworkflow.protocol.constants import (
                    STATUS_NEW,
                )

                outputStep.setStatus(STATUS_NEW)

    def _checkNewInput(self):
        """Discover new micrographs from the logical Set.

        The input Set itself is authoritative. Do not use a storage
        filename or filesystem modification time as a change detector.
        """
        micDict, self.streamClosed = self._loadInputList()

        newMics = micDict.values()
        outputStep = self._getFirstJoinStep()

        if newMics:
            dependencies = self._insertNewMicsSteps(newMics)

            if outputStep is not None:
                outputStep.addPrerequisites(*dependencies)

            self.updateSteps()

    def _loadSet(
            self,
            inputSet,
            SetClass,
            getKeyFunc,
    ):
        """Load new items from the logical Set, independently of storage."""
        self.debug(
            "Loading logical input set."
        )

        inputSet.loadAllProperties()

        newItemDict = OrderedDict()

        for item in inputSet.iterItems():
            itemKey = getKeyFunc(item)

            if itemKey not in self.micDict:
                newItemDict[itemKey] = item.clone()

        streamClosed = inputSet.isStreamClosed()

        return newItemDict, streamClosed

    def _createCtfModel(self, mic, updateSampling=False):
        """ Redefined func from the base class. """
        psd = self._getPsdPath(mic)
        ctfModel = self._ctfProgram.parseOutputAsCtf(self._getCtfOutPath(mic),
                                                     self._getCtfAvrotPath(mic),
                                                     psdFile=psd)
        ctfModel.setMicrograph(mic)
        pwutils.cleanPath(self._getCtfOutPath(mic))

        return ctfModel

    def _createOutputStep(self):
        """ Do nothing in streaming case. """
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if self.inputType == 0:
            errors.append('Movie CTF estimation is not supported yet.')

        if self.usePowerSpectra:
            mic = self._getFirstMic()
            if not hasattr(mic, "_powerSpectra"):
                errors.append("Input micrographs do not have associated power spectra.")
            else:
                self.psSampling = mic._powerSpectra.getSamplingRate()

        valueStep = round(self.stepPhaseShift.get(), 2)
        valueMin = round(self.minPhaseShift.get(), 2)
        valueMax = round(self.maxPhaseShift.get(), 2)

        if not (self.minPhaseShift < self.maxPhaseShift and
                valueStep <= (valueMax - valueMin) and
                0. <= valueMax <= 180.):
            errors.append('Wrong values for phase shift search.')

        return errors

    def _citations(self):
        return ["Mindell2003", "Rohou2015", "Elferich2024"]

    def _methods(self):
        methods = []
        if self.inputMicrographs.get() is not None:
            methods.append("We calculated the CTF of %s using CTFFind. "
                           % self.getObjectTag('inputMicrographs'))
            methods.append(self.methodsVar.get(''))
            if hasattr(self, 'outputCTF'):
                methods.append('Output CTFs: %s' % self.getObjectTag('outputCTF'))

        return methods

    # -------------------------- UTILS functions ------------------------------
    def _getMicExtra(self, mic, suffix):
        """ Return a file in extra direction with root of micFn.
        :param mic: input mic object
        :param suffix: file extension
        """
        return self._getExtraPath(pwutils.removeBaseExt(os.path.basename(
            mic.getFileName())) + '_' + suffix)

    def _getPsdPath(self, mic):
        return self._getMicExtra(mic, 'ctf.mrc')

    def _getCtfOutPath(self, mic):
        return self._getMicExtra(mic, 'ctf.txt')

    def _getCtfAvrotPath(self, mic):
        return self._getMicExtra(mic, 'ctf_avrot.txt')

    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)

        return ctf

    def _getFirstMic(self):
        """ Get first mic in the input set only once. """
        return self.getInputMicrographs().getFirstItem()
