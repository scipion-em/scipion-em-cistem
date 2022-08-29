# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.constants import BETA, SCIPION_DEBUG_NOCLEAN
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pwem.protocols import EMProtocol
from pwem.objects import CTFModel
from pwem import emlib

from .program_ctffind import ProgramCtffind
from ..convert import readCtfModelStack, parseCtffind4Output

from tomo.objects import CTFTomo, SetOfCTFTomoSeries
from tomo.protocols import ProtTsEstimateCTF


class CistemProtTsCtffind(ProtTsEstimateCTF):
    """ CTF estimation on a set of tilt series using CTFFIND4. """
    _label = 'tilt-series ctffind4'
    _devStatus = BETA
    _possibleOutputs = {'outputSetOfCTFTomoSeries': SetOfCTFTomoSeries}

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.usePowerSpectra = False
        self.useStacks = True

    # -------------------------- DEFINE param functions -----------------------
    def _initialize(self):
        ProtTsEstimateCTF._initialize(self)
        self._ctfProgram = ProgramCtffind(self)

    def _defineProcessParams(self, form):
        form.addHidden('recalculate', params.BooleanParam, default=False,
                       condition='recalculate',
                       label="Do recalculate ctf?")
        form.addHidden('continueRun', params.PointerParam, allowsNull=True,
                       condition='recalculate', label="Input previous run",
                       pointerClass='ProtTsCtffind')
        form.addHidden('sqliteFile', params.FileParam,
                       condition='recalculate',
                       allowsNull=True)
        # ctffind resamples input mics automatically
        form.addHidden('ctfDownFactor', params.FloatParam, default=1.)
        ProgramCtffind.defineProcessParams(form)

    # --------------------------- STEPS functions -----------------------------
    def processTiltSeriesStep(self, tsId):
        """ Run ctffind on a whole TS stack at once. """
        ts = self._getTiltSeries(tsId)
        tsInputFn = ts.getFirstItem().getFileName()
        tsFn = self._getTmpPath(tsId + ".mrcs")

        if pwutils.getExt(tsInputFn) in ['.mrc', '.st', '.mrcs']:
            pwutils.createAbsLink(os.path.abspath(tsInputFn), tsFn)
        else:
            ih = emlib.image.ImageHandler()
            ih.convert(tsInputFn, tsFn, emlib.DT_FLOAT)

        outputLog = self._getExtraPath(tsId + "_ctf.txt")
        outputPsd = outputLog.replace(".txt", ".mrcs")

        program, args = self._ctfProgram.getCommand(
            micFn=tsFn,
            powerSpectraPix=None,
            ctffindOut=outputLog,
            ctffindPSD=outputPsd)

        try:
            self.runJob(program, args)
            ctfResult = parseCtffind4Output(outputLog)
            ctf = CTFModel()

            for i, ti in enumerate(self._tsDict.getTiList(tsId)):
                ti.setCTF(self.getCtfTi(ctf, ctfResult, i, outputPsd))

            if not pwutils.envVarOn(SCIPION_DEBUG_NOCLEAN):
                pwutils.cleanPath(tsFn)
                pwutils.cleanPath(outputLog)

            self._tsDict.setFinished(tsId)
        except Exception as e:
            self.error(f"ERROR: Ctffind has failed for {tsFn}: {e}")

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
    def _doInsertTiltImageSteps(self):
        """ Default True, but if return False, the steps for each
        TiltImage will not be inserted. """
        return False

    def getCtfTi(self, ctf, ctfArray, tiIndex, psdStack):
        """ Parse the CTF object estimated for this Tilt-Image. """
        readCtfModelStack(ctf, ctfArray, item=tiIndex)
        ctf.setPsdFile(f"{tiIndex+1}@" + psdStack)
        ctfTomo = CTFTomo.ctfModelToCtfTomo(ctf)

        return ctfTomo
