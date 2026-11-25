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
from collections import namedtuple

from pwem.protocols import EMProtocol
from pyworkflow.object import Set
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pwem.objects import CTFModel
from pwem import emlib

from tomo.objects import CTFTomo, SetOfCTFTomoSeries, CTFTomoSeries, TiltSeries
from tomo.protocols.protocol_ts_estimate_ctf import createCtfParams

from cistem.protocols.program_ctffind import ProgramCtffind
from cistem.convert import readCtfModelStack, parseCtffindOutput
from cistem import Plugin
from cistem.protocols.protocol_ctffind import isUsableCtf

# The steps return one file here as their result file. On Continue, pyworkflow keeps a finished step
# finished only while its result files exist, so removing the ones of a tilt series that failed makes its
# steps run again.
STEP_RESULTS_DIR = 'STEP_RESULTS'
TS_STEPS = ('convert', 'process', 'output')


MRCS_EXT = ".mrcs"
# create simple, lightweight data structures similar to a class, but without the overhead of defining a full class
CistemTsCtfMd = namedtuple('CistemTsCtfMd',
                           ['ts', 'tsFn', 'outputLog', 'outputRotAvg', 'outputPsd'])


class TsCtffindOutputs(Enum):
    CTFs = SetOfCTFTomoSeries


class CistemProtTsCtffind(EMProtocol):
    """
    Estimates contrast transfer function parameters for tilt-series data using
    cisTEM Ctffind, allowing accurate characterization of image defocus,
    astigmatism, and related optical properties required for cryo-electron
    tomography reconstruction and interpretation.

    AI Generated:

    Tilt-Series Ctffind (CistemProtTsCtffind) — User Manual
        Overview

        The Tilt-Series Ctffind protocol estimates the contrast transfer
        function (CTF) of tilt-series images acquired during cryo-electron
        tomography experiments. Accurate CTF estimation is essential because
        the microscope optics modulate the transfer of structural information
        into recorded images in a resolution-dependent manner. Without proper
        correction and characterization of these optical effects, tomographic
        reconstructions may suffer from reduced contrast, distorted frequency
        information, and limited interpretability.

        In practical cryo-ET workflows, this protocol is commonly applied
        after tilt-series import and before tomographic reconstruction or
        subtomogram averaging. The protocol determines defocus values,
        astigmatism parameters, and additional diagnostic measurements for
        each tilt image, generating a consistent set of CTF estimations that
        can later be used during reconstruction, refinement, or downstream
        structural analysis.

        Biological Importance of CTF Estimation

        In cryo-electron microscopy, the CTF alters how different spatial
        frequencies are represented in recorded images. Some frequencies are
        amplified while others become attenuated or inverted. Reliable
        estimation of these effects is therefore fundamental for recovering
        high-resolution structural information.

        In tomography, the problem becomes even more important because each
        tilt image is acquired under a different viewing geometry and often
        under reduced signal-to-noise conditions. Errors in defocus estimation
        can propagate through the reconstruction process and negatively affect
        subtomogram averaging, particle localization, and interpretation of
        flexible or heterogeneous biological assemblies.

        The protocol is suitable for a wide range of biological samples,
        including membrane proteins, viral particles, organelles, cellular
        sections, and in situ tomography datasets. Accurate CTF estimation is
        especially important in high-resolution subtomogram averaging
        workflows, where small optical inaccuracies can strongly influence
        the final achievable resolution.

        Inputs and General Workflow

        The protocol accepts tilt-series datasets or previously generated CTF
        tilt-series metadata. Each tilt image within the series is analyzed
        independently while preserving the organizational structure of the
        original acquisition.

        During processing, the protocol evaluates the amplitude spectrum of
        each tilt image and determines the optical parameters that best match
        the observed signal modulation. The resulting CTF estimations are
        grouped into output CTF tilt-series objects that remain associated
        with the corresponding tilt images throughout downstream tomography
        processing.

        The workflow is designed to support both conventional tomography
        processing and more advanced cryo-ET pipelines where accurate optical
        characterization is required for refinement, dose weighting, or
        subtomogram averaging.

        Defocus and Astigmatism Estimation

        The main goal of the protocol is the estimation of objective lens
        defocus and astigmatism. In practical terms, these parameters describe
        how the microscope focus differs across image directions due to
        imperfections in the optical system.

        Biological users should interpret the estimated defocus values as a
        description of the imaging conditions rather than as intrinsic sample
        properties. Consistent defocus measurements across the tilt series
        generally indicate stable acquisition conditions, whereas large
        variations may suggest acquisition problems, specimen deformation, or
        geometric inconsistencies.

        Astigmatism measurements are also important because excessive
        astigmatism may reduce achievable resolution or introduce directional
        artifacts into tomographic reconstructions. Monitoring these values
        provides useful quality-control information during data processing.

        Tilt Geometry and Thickness Estimation

        Advanced versions of the protocol support estimation of sample tilt
        geometry and specimen thickness. These measurements can improve the
        physical realism of the CTF model and enhance parameter estimation in
        challenging tomography datasets.

        Tilt estimation attempts to determine the orientation of the specimen
        relative to the electron beam. Although computationally demanding,
        this option may improve robustness when imaging geometry significantly
        influences the observed spectra.

        Thickness estimation can provide additional physical constraints that
        improve CTF fitting quality, particularly for thicker specimens or
        cellular tomography datasets. This feature becomes especially useful
        in cryo-focused ion beam lamellae or dense in situ samples where
        multiple scattering and thickness effects become more pronounced.

        Different optimization approaches are available for thickness
        estimation. Broad systematic searches provide increased robustness at
        the expense of computation time, whereas refinement-oriented
        approaches are often faster and suitable when acquisition conditions
        are already well controlled.

        Resolution Limits and Search Parameters

        The protocol allows users to define resolution limits and defocus
        search ranges that determine the parameter space explored during CTF
        estimation. These settings strongly influence both reliability and
        computational cost.

        Lower-resolution limits help stabilize fitting by emphasizing stronger
        low-frequency signal, whereas higher-resolution limits enable more
        precise characterization when data quality permits. For noisy tilt
        images or thick samples, conservative resolution ranges are generally
        more reliable.

        Defocus search ranges should reflect the expected acquisition
        conditions. Excessively broad searches increase computation time and
        may produce unstable solutions, while overly narrow ranges can prevent
        convergence toward the correct values.

        In routine cryo-ET processing, default values are often sufficient for
        well-acquired datasets. However, difficult samples, highly tilted
        images, or low-dose conditions may require careful optimization of
        these parameters.

        Outputs and Their Interpretation

        The protocol produces a set of CTF estimations associated with the
        original tilt series. Each tilt image receives its own CTF model and
        associated diagnostic information, preserving the acquisition order
        and tomography metadata.

        Diagnostic outputs include power spectrum representations and
        rotationally averaged fitting information that can be used to evaluate
        the quality of the estimation. These diagnostics are valuable for
        identifying problematic tilt images, poor signal conditions, charging
        effects, or microscope instabilities.

        Biologically, successful CTF estimation improves the interpretability
        and consistency of tomographic reconstructions. Reliable optical
        characterization contributes directly to better subtomogram averaging,
        improved structural detail, and more accurate interpretation of
        macromolecular organization within cells and tissues.

        Practical Recommendations

        For most cryo-electron tomography datasets, it is advisable to begin
        with standard CTF estimation settings and visually inspect the
        resulting diagnostic spectra. Stable and physically reasonable
        defocus trends across tilt angles usually indicate successful
        estimation.

        For cellular tomography or thick specimens, enabling thickness-related
        refinements may improve fitting quality. However, these options also
        increase computational cost and should generally be reserved for
        datasets where higher precision is required.

        When processing highly tilted images, users should expect lower
        signal-to-noise ratios and potentially reduced fitting stability.
        Conservative resolution ranges and careful quality control become
        particularly important under these conditions.

        In high-resolution subtomogram averaging workflows, accurate CTF
        estimation is one of the most critical preprocessing steps because
        downstream refinement quality depends strongly on the reliability of
        the optical parameters assigned to each tilt image.

        Final Perspective

        Reliable CTF estimation is a foundational component of cryo-electron
        tomography data processing. By accurately characterizing microscope
        optical effects across an entire tilt series, this protocol enables
        more faithful tomographic reconstruction and more reliable biological
        interpretation. Careful parameter selection, consistent quality
        control, and awareness of specimen-specific imaging conditions are key
        elements for obtaining robust and biologically meaningful results.
    """

    _label = 'tilt-series ctffind'
    _devStatus = PROD
    _possibleOutputs = TsCtffindOutputs
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.usePowerSpectra = False
        self.useStacks = True
        self.tsCtfMdList = []
        self.inTsSet = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form, stream=False):
        form.addSection(label='Input')
        form.addParam('inputTiltSeries', params.PointerParam,
                      important=True,
                      pointerClass='SetOfTiltSeries, SetOfCTFTomoSeries',
                      label='Tilt series')
        # ctffind resamples input mics automatically
        form.addHidden('ctfDownFactor', params.FloatParam,
                       default=1.)
        ProgramCtffind.defineProcessParams(form)

        if Plugin.getActiveVersion().startswith('5'):
            form.addSection(label='Tomo')
            form.addParam('measureTilt', params.BooleanParam,
                          default=False,
                          label="Determine sample tilt?",
                          help="Measure tilt axis and angle.\n"
                               "NOTE: this will slow down computation x100 times!")

            group = form.addGroup('Thickness')
            group.addParam('measureThickness', params.BooleanParam,
                           default=False,
                           label="Determine sample thickness?\n"
                                 "NOTE: this can improve the results with a tiny "
                                 "computing time penalty!")
            group.addParam('search1D', params.BooleanParam,
                           default=True, condition='measureThickness',
                           label="Use brute force 1D search?",
                           help="CTFFIND5 will further refine thickness "
                                "and defocus by calculating the normalized "
                                "cross-correlation between the radial average "
                                "of the power spectrum (corrected for "
                                "astigmatism) and CTF, searching systematically"
                                "for the best combination of thickness in the "
                                "range of 50-400 nm in 10 nm steps, and defocus "
                                "in the range of +/-200 nm from the previously "
                                "fitted value, also in 10 nm steps.")
            group.addParam('refine2D', params.BooleanParam,
                           default=True, condition='measureThickness',
                           label="Use 2D refinement?",
                           help="CTFFIND5 will optimize thickness, defocus and "
                                "amplitude contrast using the same conjugate "
                                "gradient algorithm used in CTFFIND4 and "
                                "the normalized cross-correlation between CTF "
                                "and the 2D power spectrum as a scoring function.")
            line = group.addLine('Resolution limit for nodes (A)',
                                 condition='measureThickness')
            line.addParam('lowResNodes', params.FloatParam,
                          default=30., label='Min')
            line.addParam('highResNodes', params.FloatParam,
                          default=3., label='Max')
            group.addParam('useRoundedSquare', params.BooleanParam,
                           default=False, condition='measureThickness',
                           label="Use rounded square for nodes?")
            group.addParam('downweightNodes', params.BooleanParam,
                           default=False, condition='measureThickness',
                           label="Downweight nodes?")

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._prepareFailedForRetry()
        pIdList = []
        for mdObj in self.tsCtfMdList:
            pidConvert = self._insertFunctionStep(self.convertInputStep,
                                                  mdObj, prerequisites=[],
                                                  needsGPU=False)
            pidProcess = self._insertFunctionStep(self.processTiltSeriesStep,
                                                  mdObj, prerequisites=pidConvert,
                                                  needsGPU=False)
            pidCreateOutput = self._insertFunctionStep(self.createOutputStep,
                                                       mdObj, prerequisites=pidProcess,
                                                       needsGPU=False)
            pIdList.append(pidCreateOutput)
        self._insertFunctionStep(self.closeStep, prerequisites=pIdList, needsGPU=False)

    def _initialize(self):
        self.inTsSet = self._getInputTs()
        self._params = createCtfParams(self.inTsSet, self.windowSize.get(),
                                       self.lowRes.get(), self.highRes.get(),
                                       self.minDefocus.get(), self.maxDefocus.get())
        self._ctfProgram = ProgramCtffind(self)
        for ts in self.inTsSet.iterItems():
            outputLog = self._getExtraPath(ts.getTsId() + "_ctf.txt")
            md = CistemTsCtfMd(
                ts=ts.clone(ignoreAttrs=[]),
                tsFn=self._getTmpPath(ts.getTsId() + MRCS_EXT),
                outputLog=outputLog,
                outputRotAvg=outputLog.replace(".txt", "_avrot.txt"),
                outputPsd=outputLog.replace(".txt", MRCS_EXT)
            )
            self.tsCtfMdList.append(md)

    def convertInputStep(self, mdObj):
        tsInputFn = mdObj.ts.getFirstItem().getFileName()
        try:
            if pwutils.getExt(tsInputFn) in ['.mrc', '.st', MRCS_EXT]:
                pwutils.createAbsLink(os.path.abspath(tsInputFn), mdObj.tsFn)
            else:
                ih = emlib.image.ImageHandler()
                ih.convert(tsInputFn, mdObj.tsFn, emlib.DT_FLOAT)
        except Exception as e:
            # Not raised: the step would fail and take the whole run with it
            self.error(f"ERROR: Converting {tsInputFn} has failed: {e}")
        return self._stepResultFile('convert', mdObj)

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
        return self._stepResultFile('process', mdObj)

    def createOutputStep(self, mdObj):
        self._createOutput(mdObj)
        return self._stepResultFile('output', mdObj)

    def _readTsCtfs(self, mdObj):
        """ The CTF of every tilt image of the tilt series, or None when ctffind gave no usable result for it.

        Read in full before anything is published, so a tilt series is published whole or not at all. """
        ts = mdObj.ts
        try:
            ctfResult = parseCtffindOutput(mdObj.outputLog)
            avrotResult = parseCtffindOutput(mdObj.outputRotAvg, avrot=True)
            if ctfResult is None:
                return None
            ctf = CTFModel()
            ctfTomos = []
            for i, tiltImage in enumerate(ts.iterItems()):
                ctfTomo = self._getCtfTi(ctf, ctfResult, avrotResult, i, mdObj.outputPsd)
                ctfTomo.setIndex(tiltImage.getIndex())
                ctfTomo.setAcquisitionOrder(tiltImage.getAcquisitionOrder())
                ctfTomos.append(ctfTomo)
        except Exception as e:
            self.error(f"ERROR: The CTF of {ts.getTsId()} could not be read: {e}")
            return None
        # Some tilt images without a CTF are kept as they were; none at all is a failure
        if not any(isUsableCtf(ctfTomo) for ctfTomo in ctfTomos):
            return None
        return ctfTomos

    def _createOutput(self, mdObj):
        ts = mdObj.ts
        tsId = ts.getTsId()
        ctfTomos = self._readTsCtfs(mdObj)
        if ctfTomos is None:
            # Not published: Continue estimates it again
            self.error(f"ERROR: No CTF was estimated for {tsId}; it is left out.")
            self._failedTsIds().add(tsId)
            return

        with self._lock:
            outCtfSet = self.getOutputCtfTomoSet()
            if outCtfSet is not None and tsId in set(outCtfSet.getUniqueValues(TiltSeries.TS_ID_FIELD)):
                # Published by a previous execution of this step, e.g. before a Continue
                return

            if outCtfSet:
                outCtfSet.enableAppend()
            else:
                outCtfSet = SetOfCTFTomoSeries.create(self._getPath(),
                                                      template='ctfTomoSeries%s.sqlite')
                outCtfSet.setSetOfTiltSeries(self._getInputTs(pointer=True))
                outCtfSet.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**{self._possibleOutputs.CTFs.name: outCtfSet})
                self._defineSourceRelation(self.inTsSet, outCtfSet)

            # Generate the current CTF tomo series item
            newCTFTomoSeries = CTFTomoSeries(tsId=tsId,
                                             tiltSeriesPointer=ts)
            newCTFTomoSeries.copyInfo(ts)
            newCTFTomoSeries.setObjId(ts.getObjId())
            outCtfSet.append(newCTFTomoSeries)

            # Populate the CTF tomo series with the CTFs already read
            for ctfTomo in ctfTomos:
                newCTFTomoSeries.append(ctfTomo)

            outCtfSet.update(newCTFTomoSeries)

        self._store()

    def closeStep(self):
        # Closed in any case: an open output keeps its consumers waiting
        self._closeOutputSet()
        outCtfSet = self.getOutputCtfTomoSet()
        if self._failedTsIds() and (outCtfSet is None or len(outCtfSet) == 0):
            raise Exception("Ctffind failed for every one of the %d tilt series; see the errors above. "
                            "Continue estimates them again." % len(self._failedTsIds()))
        return self._stepResultFile('close')

    # --------------------------- Continue -----------------------------------
    def _failedTsIds(self):
        """ The tilt series of this execution with no usable CTF: in no output. """
        return self.__dict__.setdefault('_failedTsIdsSet', set())

    def _stepResultPath(self, name):
        return self._getExtraPath(STEP_RESULTS_DIR, '%s.TXT' % name)

    def _stepResultFile(self, kind, mdObj=None):
        """ The result file of a step: Continue runs the step again once it is removed. """
        pwutils.makePath(self._getExtraPath(STEP_RESULTS_DIR))
        path = self._stepResultPath('%s_%s' % (kind, mdObj.ts.getTsId()) if mdObj else kind)
        open(path, 'a').close()
        return path

    def _prepareFailedForRetry(self):
        """ Have Continue estimate again the tilt series that failed (or were never published).

        pyworkflow does not run a finished step again, so their steps, and the closing one, lose their
        result files. """
        outCtfSet = self.getOutputCtfTomoSet()
        publishedTsIds = set(outCtfSet.getUniqueValues(TiltSeries.TS_ID_FIELD)) if outCtfSet else set()
        retried = [mdObj.ts.getTsId() for mdObj in self.tsCtfMdList
                   if mdObj.ts.getTsId() not in publishedTsIds
                   and any(os.path.exists(self._stepResultPath('%s_%s' % (kind, mdObj.ts.getTsId())))
                           for kind in TS_STEPS)]
        for tsId in retried:
            pwutils.cleanPath(*[self._stepResultPath('%s_%s' % (kind, tsId)) for kind in TS_STEPS])
        if retried:
            pwutils.cleanPath(self._stepResultPath('close'))

    def allowsDelete(self, obj):
        return True

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

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

    # --------------------------- UTILS functions -----------------------------
    @staticmethod
    def _getCtfTi(ctf, ctfArray, rotAvgArray, tiIndex, psdStack):
        """ Parse the CTF object estimated for this Tilt-Image. """
        ctf, tiltAxis, tiltAngle, thickness = readCtfModelStack(ctf, ctfArray,
                                                                rotAvgArray,
                                                                item=tiIndex)
        ctf.setPsdFile(f"{tiIndex + 1}@" + psdStack)
        ctfTomo = CTFTomo.ctfModelToCtfTomo(ctf)
        if hasattr(ctf, "_rlnIceRingDensity"):
            ctfTomo.copyAttributes(ctf, '_rlnIceRingDensity')

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
