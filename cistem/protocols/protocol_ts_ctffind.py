# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *              Scipion Team (scipion@cnb.csic.es) [2]
# *
# * [1] MRC Laboratory of Molecular Biology (MRC-LMB)
# * [2] National Center of Biotechnology, CSIC, Spain
# *
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
from enum import Enum

from pwem.protocols import EMProtocol
from pyworkflow.object import Set
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pwem.objects import CTFModel
from pwem import emlib
from collections import namedtuple
from tomo.protocols.protocol_ts_estimate_ctf import createCtfParams
from .program_ctffind import ProgramCtffind
from ..convert import readCtfModelStack, parseCtffindOutput
from tomo.objects import CTFTomo, SetOfCTFTomoSeries, CTFTomoSeries

MRCS_EXT = ".mrcs"
# create simple, lightweight data structures similar to a class, but without the overhead of defining a full class
CistemTsCtfMd = namedtuple('CistemTsCtfMd', ['ts', 'tsFn', 'outputLog', 'outputPsd'])


class TsCtffindOutputs(Enum):
    CTFs = SetOfCTFTomoSeries


class CistemProtTsCtffind(EMProtocol):
    """ CTF estimation on a set of tilt series using CTFFIND. """
    _label = 'tilt-series ctffind'
    _devStatus = PROD
    _possibleOutputs = TsCtffindOutputs

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.usePowerSpectra = False
        self.useStacks = True
        self.tsCtfMdList = []
        self.inTsSet = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form, stream=False):
        form.addSection(label='Input')
        form.addParam('inputTiltSeries', params.PointerParam, important=True,
                      pointerClass='SetOfTiltSeries, SetOfCTFTomoSeries',
                      label='Input tilt series')
        form.addHidden('recalculate', params.BooleanParam,
                       default=False,
                       condition='recalculate',
                       label="Do recalculate ctf?")
        form.addHidden('continueRun', params.PointerParam,
                       allowsNull=True,
                       condition='recalculate',
                       label="Input previous run",
                       pointerClass='ProtTsCtffind')
        form.addHidden('sqliteFile', params.FileParam,
                       condition='recalculate',
                       allowsNull=True)
        # ctffind resamples input mics automatically
        form.addHidden('ctfDownFactor', params.FloatParam, default=1.)
        ProgramCtffind.defineProcessParams(form)

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self._initialize()
        pIdList = []
        for mdObj in self.tsCtfMdList:
            pidConvert = self._insertFunctionStep(self.convertInputStep, mdObj, prerequisites=[])
            pidProcess = self._insertFunctionStep(self.processTiltSeriesStep, mdObj, prerequisites=pidConvert)
            pidCreateOutput = self._insertFunctionStep(self.createOutputStep, mdObj, prerequisites=pidProcess)
            pIdList.append(pidCreateOutput)
        self._insertFunctionStep(self.closeStep, prerequisites=pIdList)

    def _initialize(self):
        self.inTsSet = self._getInputTs()
        self._params = createCtfParams(self.inTsSet, self.windowSize.get(), self.lowRes.get(), self.highRes.get(),
                                       self.minDefocus.get(), self.maxDefocus.get())
        self._ctfProgram = ProgramCtffind(self)
        for ts in self.inTsSet.iterItems():
            outputLog = self._getExtraPath(ts.getTsId() + "_ctf.txt")
            md = CistemTsCtfMd(
                ts=ts.clone(ignoreAttrs=[]),
                tsFn=self._getTmpPath(ts.getTsId() + MRCS_EXT),
                outputLog=outputLog,
                outputPsd=outputLog.replace(".txt", MRCS_EXT)
            )
            self.tsCtfMdList.append(md)

    @staticmethod
    def convertInputStep(mdObj):
        tsInputFn = mdObj.ts.getFirstItem().getFileName()
        if pwutils.getExt(tsInputFn) in ['.mrc', '.st', MRCS_EXT]:
            pwutils.createAbsLink(os.path.abspath(tsInputFn), mdObj.tsFn)
        else:
            ih = emlib.image.ImageHandler()
            ih.convert(tsInputFn, mdObj.tsFn, emlib.DT_FLOAT)

    def processTiltSeriesStep(self, mdObj):
        """ Run ctffind on a whole TS stack at once. """
        program, args = self._ctfProgram.getCommand(
            micFn=mdObj.tsFn,
            powerSpectraPix=None,
            ctffindOut=mdObj.outputLog,
            ctffindPSD=mdObj.outputPsd)

        try:
            self.runJob(program, args)
        except Exception as e:
            self.error(f"ERROR: Ctffind has failed for {mdObj.tsFn}: {e}")

    def createOutputStep(self, mdObj):
        outCtfSet = self.getOutputCtfTomoSet()

        if outCtfSet:
            outCtfSet.enableAppend()
        else:
            outCtfSet = SetOfCTFTomoSeries.create(self._getPath(), template='ctfTomoSeriess%s.sqlite')
            outCtfSet.setSetOfTiltSeries(self.inTsSet)
            outCtfSet.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{self._possibleOutputs.CTFs.name: outCtfSet})
            self._defineSourceRelation(self.inTsSet, outCtfSet)

        ts = mdObj.ts
        outputLog = mdObj.outputLog

        # Generate the current CTF tomo series item
        newCTFTomoSeries = CTFTomoSeries()
        newCTFTomoSeries.copyInfo(ts)
        newCTFTomoSeries.setTiltSeries(ts)
        newCTFTomoSeries.setObjId(ts.getObjId())
        newCTFTomoSeries.setTsId(ts.getTsId())
        outCtfSet.append(newCTFTomoSeries)

        # Generate the ti CTF and populate the corresponding CTF tomo series
        ctfResult = parseCtffindOutput(outputLog)
        ctf = CTFModel()
        for i, tiltImage in enumerate(ts.iterItems()):
            ctfTomo = self._getCtfTi(ctf, ctfResult, i, mdObj.outputPsd)
            ctfTomo.setIndex(tiltImage.getIndex())
            ctfTomo.setAcquisitionOrder(tiltImage.getAcquisitionOrder())
            tiltImage.setCTF(ctfTomo)
            newCTFTomoSeries.append(ctfTomo)

        outCtfSet.update(newCTFTomoSeries)
        self._store()

    def closeStep(self):
        outCtfSet = self.getOutputCtfTomoSet()
        outCtfSet.setStreamState(Set.STREAM_CLOSED)
        outCtfSet.write()
        self._store()

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if self.lowRes.get() > 50:
            errors.append("Minimum resolution cannot be > 50A.")

        valueStep = round(self.stepPhaseShift.get(), 2)
        valueMin = round(self.minPhaseShift.get(), 2)
        valueMax = round(self.maxPhaseShift.get(), 2)

        if not (self.minPhaseShift < self.maxPhaseShift and
                valueStep <= (valueMax - valueMin) and
                0. <= valueMax <= 180.):
            errors.append('Wrong values for phase shift search.')

        return errors

    def _citations(self):
        return ["Mindell2003", "Rohou2015"]

    # --------------------------- UTILS functions -----------------------------
    @staticmethod
    def _getCtfTi(ctf, ctfArray, tiIndex, psdStack):
        """ Parse the CTF object estimated for this Tilt-Image. """
        readCtfModelStack(ctf, ctfArray, item=tiIndex)
        ctf.setPsdFile(f"{tiIndex + 1}@" + psdStack)
        ctfTomo = CTFTomo.ctfModelToCtfTomo(ctf)
        return ctfTomo

    def _getInputTs(self, pointer=False):
        if isinstance(self.inputTiltSeries.get(), SetOfCTFTomoSeries):
            return self.inputTiltSeries.get().getSetOfTiltSeries(pointer=pointer)
        return self.inputTiltSeries.get() if not pointer else self.inputTiltSeries

    def getOutputCtfTomoSet(self):
        return getattr(self, TsCtffindOutputs.CTFs.name, None)

    def getCtfParamsDict(self):
        """ Return a copy of the global params dict,
        to avoid overwriting values. """
        return self._params


