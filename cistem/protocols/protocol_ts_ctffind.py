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
from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pwem.protocols import EMProtocol

from .program_ctffind import ProgramCtffind

from tomo.objects import CTFTomo
from tomo.protocols import ProtTsEstimateCTF


class CistemProtTsCtffind(ProtTsEstimateCTF):
    """ CTF estimation on a set of tilt series using CTFFIND4. """
    _label = 'tilt-series ctffind4'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.usePowerSpectra = False

    # -------------------------- DEFINE param functions -----------------------
    def _initialize(self):
        ProtTsEstimateCTF._initialize(self)
        self._ctfProgram = ProgramCtffind(self)

    def _defineProcessParams(self, form):
        form.addParam('recalculate', params.BooleanParam, default=False,
                      condition='recalculate',
                      label="Do recalculate ctf?")
        form.addParam('continueRun', params.PointerParam, allowsNull=True,
                      condition='recalculate', label="Input previous run",
                      pointerClass='ProtTsCtffind')
        form.addHidden('sqliteFile', params.FileParam,
                       condition='recalculate',
                       allowsNull=True)
        # ctffind resamples input mics automatically
        form.addHidden('ctfDownFactor', params.FloatParam, default=1.)
        ProgramCtffind.defineProcessParams(form)

    # --------------------------- STEPS functions -----------------------------
    def _estimateCtf(self, workingDir, tiFn, ti):
        try:
            tsId = ti.getTsId()
            pwutils.makePath(self._getExtraPath(tsId))

            outputLog = os.path.join(workingDir, 'output-log.txt')
            outputPsd = os.path.join(workingDir, self.getPsdName(ti))

            program, args = self._ctfProgram.getCommand(
                micFn=tiFn,
                powerSpectraPix=None,
                ctffindOut=outputLog,
                ctffindPSD=outputPsd)

            self.runJob(program, args)

            # Move files we want to keep
            pwutils.moveFile(outputPsd, self._getExtraPath(tsId))
            pwutils.moveFile(outputPsd.replace('.mrc', '.txt'),
                             self._getTmpPath())
            pwutils.moveFile(outputPsd.replace('.mrc', '_avrot.txt'),
                             self._getExtraPath(tsId))
        except:
            print("ERROR: Ctffind has failed for %s" % tiFn)

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
    def getPsdName(self, ti):
        return '%s_PSD.mrc' % self.getTiPrefix(ti)

    def getCtf(self, ti):
        """ Parse the CTF object estimated for this Tilt-Image. """
        psd = self.getPsdName(ti)
        outCtf = self._getTmpPath(psd.replace('.mrc', '.txt'))
        ctfModel = self._ctfProgram.parseOutputAsCtf(outCtf,
                                                     psdFile=self._getExtraPath(ti.getTsId(), psd))
        ctfTomo = CTFTomo.ctfModelToCtfTomo(ctfModel)

        return ctfTomo
