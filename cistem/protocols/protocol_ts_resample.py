# **************************************************************************
# *
# * Authors:     Ricardo D. Righetto (ricardo.righetto@unibas.ch)
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * University of Basel
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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

# Inspired by protocol_resizeTS.py from scipion-em-xmipptomo:
# https://github.com/I2PC/scipion-em-xmipptomo/blob/devel/xmipptomo/protocols/protocol_resizeTS.py

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, IntParam, Positive
import pyworkflow.utils as pwutils
from pwem.protocols import EMProtocol
from tomo.objects import SetOfTiltSeries
from tomo.protocols import ProtTomoBase

from cistem import Plugin

OUTPUT_TS_NAME = 'resampledTiltSeries'


class CistemProtTsResample(EMProtocol, ProtTomoBase):
    # the link to cistem.org warns about a security problem
    """
    Resamples cryo-electron tomography tilt series using Fourier-space
    cropping or padding in order to change image dimensions while
    minimizing aliasing artifacts. More info: https://cistem.org

    AI Generated:

    Tilt Series Resampling (CistemProtTsResample) — User Manual
        Overview

        The Tilt Series Resampling protocol changes the image dimensions
        of cryo-electron tomography tilt series using Fourier-based
        interpolation methods provided by cisTEM. Its main purpose is to
        resize tilt images while preserving structural information and
        avoiding the aliasing artifacts commonly associated with simpler
        interpolation or pixel averaging approaches.

        In cryo-electron tomography workflows, resampling is frequently
        performed to adapt datasets to computational requirements,
        reconstruction strategies, or downstream processing tools.
        Biological users often employ this operation to reduce dataset
        size for rapid exploratory analysis or to standardize pixel sizes
        across multiple datasets before comparative studies or subtomogram
        averaging.

        Unlike conventional spatial-domain binning, Fourier resampling
        preserves frequency information more accurately. This makes the
        protocol especially useful in workflows where maintaining signal
        fidelity is important for downstream interpretation or alignment
        stability.

        Inputs and General Workflow

        The protocol requires a set of aligned or unaligned tilt series
        as input. Each tilt image is resized to new user-defined X and Y
        dimensions while preserving the overall geometry and structure of
        the dataset. The resulting output is a new tilt series set with
        updated sampling properties consistent with the applied rescaling.

        From a biological perspective, reducing image dimensions can
        substantially accelerate computationally intensive tomography
        operations such as tilt-series alignment, tomogram reconstruction,
        denoising, or subtomogram averaging. This is particularly useful
        during exploratory processing stages where rapid feedback is more
        important than preserving the maximum achievable resolution.

        Conversely, increasing the image size through Fourier padding may
        facilitate visualization or compatibility with downstream tools,
        although it does not introduce new structural information beyond
        the original experimental signal.

        Fourier-Based Resampling Strategy

        The protocol uses Fourier cropping or Fourier padding to modify
        image dimensions. Fourier cropping effectively removes the highest
        spatial frequencies, producing a downsampled representation that
        retains the most biologically relevant large-scale structural
        information while reducing noise and computational burden.

        Fourier padding enlarges the image grid without altering the
        intrinsic resolution content of the data. This operation can be
        useful for visualization purposes, coordinate consistency, or
        compatibility with reconstruction algorithms requiring specific
        image dimensions.

        Compared with direct interpolation methods, Fourier resampling
        provides smoother and more physically consistent results because
        it operates directly in frequency space. In practical cryo-ET
        workflows, this generally leads to improved preservation of
        structural continuity and reduced interpolation artifacts.

        Sampling Rate Considerations

        Resampling changes the effective sampling rate of the tilt series.
        When images are downsampled, the apparent pixel size increases
        proportionally. This reduces the maximum recoverable resolution
        but often improves signal-to-noise ratio and computational
        efficiency. Biological users should carefully balance these
        factors according to the objectives of the project.

        For rapid screening, particle localization, or low-resolution
        structural interpretation, aggressive downsampling is often
        acceptable and highly efficient. However, for high-resolution
        subtomogram averaging or detailed membrane analysis, preserving
        sufficient sampling density remains critical.

        Upsampling through Fourier padding does not improve intrinsic
        resolution and should not be interpreted as adding new biological
        detail. Instead, it primarily affects the numerical representation
        and grid dimensions of the dataset.

        Outputs and Their Interpretation

        The protocol produces a new set of tilt series with updated image
        dimensions and appropriately adjusted sampling information. The
        geometric consistency of the tilt series is preserved, making the
        outputs suitable for subsequent tomography workflows including
        alignment, reconstruction, segmentation, and subtomogram analysis.

        Because Fourier-based resampling preserves low-frequency signal
        effectively, the resulting datasets generally maintain the major
        structural features required for biological interpretation even
        after significant size reduction.

        Practical Recommendations

        In routine cryo-ET processing, moderate downsampling is often an
        excellent strategy for accelerating exploratory analyses and
        reducing storage requirements. Users commonly perform initial
        alignment and reconstruction at reduced sampling before returning
        to full-resolution datasets for final refinement.

        Care should be taken not to downsample excessively when studying
        small macromolecular complexes, fine membrane features, or subtle
        conformational variability. Excessive reduction in sampling may
        remove biologically meaningful high-frequency information that
        cannot be recovered later.

        When combining datasets from different microscopes or acquisition
        conditions, resampling can also help standardize pixel sizes and
        simplify comparative analyses across experiments.

        Final Perspective

        In cryo-electron tomography, image resampling is not merely a
        technical preprocessing step but an important compromise between
        computational efficiency, storage demands, and preservation of
        biological information. Fourier-based resizing provides a robust
        and artifact-minimizing strategy that supports both exploratory
        and high-quality tomography workflows while maintaining reliable
        structural consistency across the dataset.
    """

    _label = 'resample tilt series'
    _possibleOutputs = {OUTPUT_TS_NAME: SetOfTiltSeries}
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        #     Example 'resample' input:
        # 
        # resample
        #  
        #  
        #         **   Welcome to Resample   **
        #  
        #              Version : 1.00
        #             Compiled : Dec  2 2017
        #                 Mode : Interactive
        #  
        # Input image file name [ts2_L1G1-dose_filt.rec]   : ts2_L1G1-dose_filt.st
        # Output image file name [ts2_L1G1_72_d3-bin8.rec] : ts2_L1G1-dose_filt-bin4.st
        # Is the input a volume [YES]                        : NO
        # New X-Size [464]                                   : 928
        # New Y-Size [464]                                   : 928

        # You need a params to belong to a section:
        form.addSection(label=pwutils.Message.LABEL_INPUT)
        form.addParam('inputSetOfTiltSeries', PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt series')

        form.addParam('newXsize', IntParam,
                      default=512, validators=[Positive],
                      label='New X-Size (px)',
                      help='Images will be rescaled to this size in X dimension (pixels)')
        form.addParam('newYsize', IntParam,
                      default=512, validators=[Positive],
                      label='New Y-Size (px)',
                      help='Images will be rescaled to this size in Y dimension (pixels)')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():

            self._insertFunctionStep(self.runTsResample,
                                     ts.getFirstItem().getFileName(),
                                     needsGPU=False)

        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def runTsResample(self, tsFile: str):
        prog = Plugin.getProgram('resample')

        paramDict = {
            'tsFile': tsFile,
            'tsOutName': self._getOutputFn(tsFile),
            'newXsize': self.newXsize,
            'newYsize': self.newYsize,
        }

        # Arguments to the resample command defined in the plugin initialization:
        args = """   << eof
%(tsFile)s
%(tsOutName)s
NO
%(newXsize)d
%(newYsize)d
eof\n
"""
        self.runJob(prog, args % paramDict)

    def createOutputStep(self):
        inputTs = self.inputSetOfTiltSeries.get()
        outTsSet = self._createSetOfTiltSeries()
        outTsSet.copyInfo(inputTs)
        outTsSet.setSamplingRate(self._getOutputSampling())
        outTsSet.copyItems(inputTs, updateTiCallback=self._updateTi)

        self._defineOutputs(**{OUTPUT_TS_NAME: outTsSet})
        self._defineTransformRelation(self.inputSetOfTiltSeries, outTsSet)

    # --------------------------- INFO functions -----------------------------------
    def _citations(self):
        return ['Grant2018']

    # --------------------------- UTILS functions -------------------------------
    def _updateTi(self, j, ts, ti, tsOut, tiOut):
        fn = ti.getFileName()
        tiOut.setFileName(self._getOutputFn(fn))

    def _getOutputFn(self, tomoFile):
        tomoBaseName = pwutils.removeBaseExt(tomoFile)
        tomoExt = pwutils.getExt(tomoFile)
        output = self._getExtraPath(tomoBaseName + '_resampled' + tomoExt)

        return output

    def _getOutputSampling(self):
        ts = self.inputSetOfTiltSeries.get()
        oldSamplingRate = ts.getSamplingRate()
        oldXsize = ts.getFirstItem().getXDim()

        return oldSamplingRate * oldXsize / self.newXsize.get()
