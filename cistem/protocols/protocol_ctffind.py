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

import pyworkflow as pw

from cistem.constants import *
from .program_ctffind import ProgramCtffind


class CistemProtCTFFind(pw.em.ProtCTFMicrographs):
    """
    Estimates CTF for a set of micrographs/movies with ctffind4.
    
    To find more information about ctffind4 go to:
    http://ctffind.janelia.org/ctffind4
    """
    _label = 'ctffind4'

    def _defineParams(self, form):
        ProgramCtffind.defineFormParams(form)
        self._defineStreamingParams(form)

    def _defineCtfParamsDict(self):
        pw.em.ProtCTFMicrographs._defineCtfParamsDict(self)
        self._ctfProgram = ProgramCtffind(self)

    # -------------------------- STEPS functions ------------------------------
    def _doCtfEstimation(self, mic, **kwargs):
        """ Run ctffind with required parameters """
        micFn = mic.getFileName()
        micDir = self._getTmpPath('mic_%06d' % mic.getObjId())
        # Create micrograph dir
        pw.utils.makePath(micDir)
        micFnMrc = os.path.join(micDir, pw.utils.replaceBaseExt(micFn, 'mrc'))

        ih = pw.em.ImageHandler()

        if not ih.existsLocation(micFn):
            raise Exception("Missing input micrograph: %s" % micFn)

        ih.convert(micFn, micFnMrc, pw.em.DT_FLOAT)

        try:
            program, args = self._ctfProgram.getCommand(
                micFn=micFnMrc,
                ctffindOut=self._getCtfOutPath(mic),
                ctffindPSD=self._getPsdPath(mic),
                **kwargs)
            self.runJob(program, args)

            pw.utils.cleanPath(micDir)

        except:
            print("ERROR: Ctffind has failed for %s" % micFnMrc)
            import traceback
            traceback.print_exc()

    def _estimateCTF(self, mic, *args):
        self._doCtfEstimation(mic)

    def _reEstimateCTF(self, mic, ctf):
        self._doCtfEstimation(mic, **self._getRecalCtfParamsDict(ctf))

    def _createCtfModel(self, mic, updateSampling=False):
        psd = self._getPsdPath(mic)
        ctfModel = self._ctfProgram.parseOutputAsCtf(self._getCtfOutPath(mic),
                                                     psdFile=psd)
        ctfModel.setMicrograph(mic)

        return ctfModel

    def _createOutputStep(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if self.inputType == 0:
            errors.append('Movie CTF estimation is not ready yet.')

        valueStep = round(self.stepPhaseShift.get(), 2)
        valueMin = round(self.minPhaseShift.get(), 2)
        valueMax = round(self.maxPhaseShift.get(), 2)

        if not (self.minPhaseShift < self.maxPhaseShift and
                valueStep <= (valueMax-valueMin) and
                0. <= valueMax <= 180.):
            errors.append('Wrong values for phase shift search.')

        return errors

    def _citations(self):
        return ["Mindell2003", "Rohou2015"]

    def _methods(self):
        if self.inputMicrographs.get() is None:
            return ['Input micrographs not available yet.']
        methods = ("We calculated the CTF of %s using CTFFind. "
                   % self.getObjectTag('inputMicrographs'))
        methods += self.methodsVar.get('')
        methods += 'Output CTFs: %s' % self.getObjectTag('outputCTF')

        return [methods]

    # -------------------------- UTILS functions ------------------------------
    def _getRecalCtfParamsDict(self, ctfModel):
        values = map(float, ctfModel.getObjComment().split())
        sampling = self.inputMicrographs.get().getSamplingRate()
        return {
            'step_focus': 500.0,
            'lowRes': sampling / values[3],
            'highRes': sampling / values[4],
            'minDefocus': min([values[0], values[1]]),
            'maxDefocus': max([values[0], values[1]])
        }

    def _getMicExtra(self, mic, suffix):
        """ Return a file in extra direction with root of micFn. """
        return self._getExtraPath(pw.utils.removeBaseExt(os.path.basename(
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
        ctf = pw.em.CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)

        return ctf

    def _summary(self):
        summary = pw.em.ProtCTFMicrographs._summary(self)
        return summary
