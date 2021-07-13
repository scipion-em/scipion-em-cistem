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

import os

import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD
from pwem.protocols import ProtCTFMicrographs
from pwem.objects import CTFModel
from pwem import emlib

from .program_ctffind import ProgramCtffind


class CistemProtCTFFind(ProtCTFMicrographs):
    """ Estimate CTF for a set of micrographs with ctffind4.
    
    To find more information about ctffind4 visit:
    https://grigoriefflab.umassmed.edu/ctffind4
    """
    _label = 'ctffind4'
    _devStatus = PROD

    def _defineParams(self, form):
        ProgramCtffind.defineInputParams(form)
        ProgramCtffind.defineProcessParams(form)
        self._defineStreamingParams(form)

    def _defineCtfParamsDict(self):
        ProtCTFMicrographs._defineCtfParamsDict(self)
        self._ctfProgram = ProgramCtffind(self)

    # -------------------------- STEPS functions ------------------------------
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
        # Create micrograph dir
        pwutils.makePath(micDir)
        micFnMrc = os.path.join(micDir, pwutils.replaceBaseExt(micFn, 'mrc'))

        ih = emlib.image.ImageHandler()

        if not os.path.exists(micFn):
            raise FileNotFoundError("Missing input micrograph: %s" % micFn)

        if micFn.endswith('.mrc'):
            pwutils.createAbsLink(os.path.abspath(micFn), micFnMrc)
        else:
            ih.convert(micFn, micFnMrc, emlib.DT_FLOAT)

        try:
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
        with open(file, "r") as fh:
            for line in fh.readlines():
                if line.startswith("Error"):
                    return line.replace("Error:", "")
        return e

    def _estimateCTF(self, mic, *args):
        """ Redefined func from the base class. """
        self._doCtfEstimation(mic)

    def _reEstimateCTF(self, mic, ctf):
        """ Redefined func from the base class. """
        self._doCtfEstimation(mic, **self._getRecalCtfParamsDict(ctf))

    def _createCtfModel(self, mic, updateSampling=False):
        """ Redefined func from the base class. """
        psd = self._getPsdPath(mic)
        ctfModel = self._ctfProgram.parseOutputAsCtf(self._getCtfOutPath(mic),
                                                     psdFile=psd)
        ctfModel.setMicrograph(mic)

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

    def _methods(self):
        methods = []
        if self.inputMicrographs.get() is not None:
            methods.append("We calculated the CTF of %s using CTFFind. "
                           % self.getObjectTag('inputMicrographs'))
            methods.append(self.methodsVar.get(''))
            if getattr(self, 'outputCTF'):
                methods.append('Output CTFs: %s' % self.getObjectTag('outputCTF'))

        return methods

    # -------------------------- UTILS functions ------------------------------
    def _getRecalCtfParamsDict(self, ctfModel):
        """ Update values from user-adjusted params from the Java GUI.
        :param ctfModel: input CTF model object
        """
        values = [float(x) for x in ctfModel.getObjComment().split()]
        sampling = self.inputMicrographs.get().getSamplingRate()
        return {
            'step_focus': 500.0,
            'lowRes': sampling / values[3],
            'highRes': sampling / values[4],
            'minDefocus': min([values[0], values[1]]),
            'maxDefocus': max([values[0], values[1]])
        }

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

    def _parseOutput(self, filename):
        """ Try to find the output estimation parameters
        from filename. It searches for a first line without #.
        """
        return self._ctfProgram.parseOutput(filename)

    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)

        return ctf

    def _getFirstMic(self):
        """ Get first mic in the input set only once. """
        return self.getInputMicrographs().getFirstItem()
