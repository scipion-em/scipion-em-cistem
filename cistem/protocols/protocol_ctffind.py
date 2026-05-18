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
from pyworkflow.object import Boolean
from pwem.protocols import ProtCTFMicrographs
from pwem.objects import CTFModel
from pwem import emlib

from .program_ctffind import ProgramCtffind


class CistemProtCTFFind(ProtCTFMicrographs):
    """
    Estimates the Contrast Transfer Function (CTF) parameters of cryo-EM
    micrographs using CTFFIND. The protocol evaluates the optical effects
    introduced by the electron microscope and determines parameters such as
    defocus, astigmatism, and phase shift, which are essential for accurate
    downstream reconstruction and image processing. More info:
        https://grigoriefflab.umassmed.edu/ctffind4

    AI Generated:

    CTFFIND CTF Estimation (CistemProtCTFFind) — User Manual
        Overview

        The CTFFIND CTF Estimation protocol determines the optical transfer
        characteristics of cryo-EM micrographs by analyzing their frequency
        patterns. In cryo-electron microscopy, the Contrast Transfer Function
        describes how different spatial frequencies are affected by the imaging
        conditions of the microscope. Accurate estimation of these parameters is
        one of the most important preprocessing steps because nearly all later
        stages of reconstruction depend on reliable CTF correction.

        For biological users, this protocol is typically applied immediately
        after motion correction and before particle picking or classification.
        The resulting defocus and astigmatism measurements are critical for
        achieving high-resolution reconstructions and for evaluating the overall
        quality of the dataset.

        Inputs and General Workflow

        The protocol operates on a collection of micrographs or associated power
        spectra derived from them. In most workflows, users provide motion-
        corrected micrographs directly, allowing the protocol to compute the
        corresponding frequency-domain information internally. Alternatively,
        externally generated power spectra may be supplied when preprocessing
        pipelines already include spectral estimation steps.

        During execution, the protocol evaluates the oscillatory signal present
        in the Fourier transform of each micrograph and identifies the optical
        parameters that best explain the observed Thon ring patterns. The
        resulting CTF models are then associated with the corresponding
        micrographs and become available for downstream cryo-EM processing.

        Biological Importance of Accurate CTF Estimation

        Correct CTF estimation directly influences the interpretability and
        resolution of cryo-EM reconstructions. Poorly estimated defocus values
        propagate errors throughout the processing workflow, reducing map
        quality and potentially introducing artifacts into reconstructed
        structures.

        In practical biological studies, datasets often contain micrographs
        acquired under slightly different optical conditions. The protocol
        allows each image to be analyzed independently so that local variations
        in defocus and astigmatism can be accurately captured. This becomes
        particularly important in high-resolution single-particle analysis,
        where small inaccuracies can significantly limit the final resolution.

        Phase Shift Estimation

        The protocol supports phase shift estimation for datasets collected
        using phase plates. In these experiments, the microscope intentionally
        modifies image contrast to improve visualization of weakly scattering
        biological specimens. Estimating the phase shift correctly is essential
        because inaccurate values may distort the recovered structural signal.

        Biological users should define realistic phase shift search ranges that
        reflect the acquisition conditions used during data collection. Very
        broad or inconsistent search ranges may reduce stability and lead to
        unreliable solutions.

        Use of Power Spectra

        In some workflows, users may choose to estimate CTF parameters directly
        from precomputed power spectra instead of raw micrographs. This approach
        can be useful in automated acquisition systems or specialized facility
        pipelines where spectral preprocessing has already been optimized.

        However, the quality of the provided spectra strongly determines the
        robustness of the estimation. Poor spectral normalization, excessive
        masking, or contamination may interfere with accurate fitting of the
        Thon rings.

        Interpretation of Results

        After execution, the protocol generates a CTF model for each processed
        micrograph. These models contain the estimated optical parameters and
        associated diagnostic information used throughout the remainder of the
        cryo-EM workflow.

        From a biological perspective, users should inspect the consistency of
        the estimated defocus values across the dataset. Large variations may
        indicate acquisition instability, ice thickness heterogeneity, charging
        effects, or contamination. Likewise, unusually high astigmatism values
        can reflect microscope misalignment or data collection problems.

        The generated power spectrum diagnostic images are particularly useful
        for visual quality control. Well-defined and continuous Thon rings
        generally indicate reliable estimation and good micrograph quality,
        whereas weak or discontinuous rings may reveal drift, poor vitrification,
        contamination, or low signal-to-noise conditions.

        Streaming and High-Throughput Workflows

        The protocol is compatible with streaming-oriented cryo-EM workflows,
        allowing CTF estimation to proceed as new micrographs become available.
        This capability is especially important in modern automated facilities,
        where rapid feedback during acquisition helps users detect imaging
        problems early and optimize microscope conditions in real time.

        In high-throughput environments, streaming estimation can significantly
        accelerate decision-making by enabling immediate assessment of defocus
        distributions, phase shift stability, and overall dataset quality while
        data collection is still ongoing.

        Practical Recommendations

        In routine biological practice, it is generally advisable to begin with
        conservative default parameters and carefully inspect the resulting
        diagnostic spectra before processing the entire dataset. Consistent Thon
        ring visibility across multiple micrographs is often the best indicator
        of stable acquisition conditions.

        For phase-plate datasets, users should pay particular attention to the
        selected phase shift search limits and verify that estimated values are
        biologically and experimentally reasonable. Excessively noisy
        micrographs or images containing crystalline ice contamination may
        produce unstable estimations and are often best excluded from further
        processing.

        When working at very high resolution, accurate calibration of pixel size,
        voltage, spherical aberration, and amplitude contrast becomes especially
        important because small parameter inaccuracies can affect downstream
        refinement quality.

        Final Perspective

        CTF estimation is one of the foundational steps of cryo-EM image
        processing because it defines how microscope optics have altered the
        recorded biological signal. Reliable estimation enables accurate
        correction of imaging artifacts and establishes the basis for high-
        resolution structural interpretation.

        For most cryo-EM projects, careful inspection of CTF quality and
        thoughtful interpretation of the resulting parameters are as important
        as the numerical estimation itself. Stable optical measurements,
        consistent Thon ring patterns, and biologically plausible defocus values
        are key indicators of a high-quality dataset suitable for reliable
        downstream reconstruction.
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
                if "Error:" in line:
                    return line.split("Error:")[-1]
        return e

    def _estimateCTF(self, mic, *args):
        """ Redefined func from the base class. """
        self._doCtfEstimation(mic)

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
