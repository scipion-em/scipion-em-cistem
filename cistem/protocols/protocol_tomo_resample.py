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
import os.path

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, IntParam, Positive
import pyworkflow.utils as pwutils
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfTomograms

from cistem import Plugin

OUTPUT_TOMO_NAME = 'resampledTomos'


class CistemProtTomoResample(EMProtocol, ProtTomoBase):
    #the link to cistem.org warns about a security problem
    """
    Resamples tomograms using Fourier-based cropping and padding methods
    implemented in cisTEM. The protocol allows tomographic volumes to be
    resized while preserving structural fidelity and minimizing aliasing
    artifacts that are commonly introduced by conventional interpolation or
    binning approaches. This operation is particularly useful in cryo-electron
    tomography workflows where tomograms must be standardized to compatible
    voxel sizes or dimensions before downstream processing, visualization,
    segmentation, subtomogram averaging, or data sharing. More info:
    https://cistem.org

    AI Generated:

    Resample Tomogram (CistemProtTomoResample) - User Manual
        Overview

        The Resample Tomogram protocol changes the dimensions and effective
        sampling of tomographic reconstructions while preserving the underlying
        frequency information as accurately as possible. Instead of relying on
        direct interpolation in real space, the protocol performs Fourier-space
        resizing, which generally produces cleaner and more biologically reliable
        results. This makes the protocol especially valuable in cryo-electron
        tomography workflows where tomograms are frequently adapted for different
        computational or visualization requirements.

        In practical biological research, tomograms are often generated at very
        large sizes that can become difficult to process efficiently. Resampling
        allows users to reduce computational cost, improve visualization speed,
        or standardize datasets acquired under different imaging conditions.
        Conversely, users may also increase dimensions when compatibility with
        downstream software requires a specific voxel grid.

        Biological Context and Typical Applications

        Tomogram resampling is commonly performed before segmentation, template
        matching, subtomogram extraction, neural-network analysis, or visualization
        in external software packages. Large tomograms may exceed available GPU or
        memory resources, making reduced dimensions necessary for practical
        processing. In these situations, controlled Fourier resampling helps retain
        structural interpretability while reducing data size.

        Another common use case involves harmonizing datasets acquired at different
        pixel sizes. When combining tomograms from multiple microscopes, acquisition
        sessions, or processing pipelines, maintaining consistent sampling becomes
        important for quantitative comparison and downstream averaging procedures.

        From a biological perspective, users should remember that resampling changes
        the effective voxel spacing of the reconstructed specimen. Any interpretation
        involving distances, membrane thickness, molecular dimensions, or spatial
        organization should therefore always consider the updated sampling rate after
        processing.

        Choosing the Output Dimensions

        The protocol allows independent specification of the X, Y, and Z dimensions
        of the output tomogram. This flexibility is useful because tomographic
        datasets are often anisotropic, particularly along the Z direction where the
        missing wedge and acquisition geometry influence reconstruction quality.

        In many workflows, reducing the X and Y dimensions is sufficient to decrease
        storage and computational burden. Z resampling may also be beneficial when
        preparing tomograms for machine-learning pipelines or visualization software
        that performs better with isotropic or standardized volumes.

        Users should select dimensions that remain biologically meaningful for the
        structures under investigation. Excessive downsampling can erase fine
        structural details such as membranes, filaments, ribosomes, or macromolecular
        assemblies. Moderate reductions are generally safer and often provide a good
        balance between computational efficiency and structural preservation.

        Fourier-Based Resampling and Image Quality

        The protocol uses Fourier-space resizing to minimize aliasing artifacts and
        preserve frequency information more faithfully than simple interpolation
        methods. This distinction is biologically important because interpolation-based
        resizing can introduce artificial smoothing or distortions that compromise
        interpretation of subtle structural features.

        In cryo-electron tomography, maintaining reliable spatial frequencies is
        particularly important for downstream analyses involving correlation,
        classification, segmentation, or subtomogram averaging. Fourier-based
        resampling therefore provides a safer strategy for preparing tomograms while
        reducing the risk of introducing processing artifacts.

        Nevertheless, users should recognize that any reduction in sampling decreases
        the maximum recoverable spatial detail. Even when performed carefully, strong
        downsampling inevitably limits the interpretability of high-resolution features.

        Outputs and Interpretation

        The protocol produces a new set of tomograms with updated dimensions and
        corresponding sampling information. The resulting datasets preserve the
        identity and metadata of the original tomograms while adapting them to the
        requested output geometry.

        Because the voxel size changes after resampling, downstream measurements and
        visualization settings should always use the updated sampling information.
        Failure to account for the new voxel spacing may lead to incorrect biological
        interpretation or inaccurate dimensional measurements.

        Practical Recommendations

        For most biological workflows, moderate downsampling is often sufficient to
        accelerate processing without substantially affecting interpretability. This
        is particularly useful during exploratory analysis, manual segmentation, or
        rapid visualization tasks.

        When preparing tomograms for quantitative structural analysis, users should
        avoid aggressive reductions that may erase biologically meaningful detail.
        If high-resolution interpretation is required, retaining the original sampling
        or using only mild resampling is generally preferable.

        It is also advisable to maintain consistency across related datasets. Using
        standardized voxel sizes simplifies comparison between tomograms and improves
        compatibility with downstream computational workflows.

        Final Perspective

        Tomogram resampling is not simply a technical resizing operation but an
        important preparation step that influences visualization quality, computational
        performance, and biological interpretation. Careful selection of output
        dimensions and thoughtful consideration of the required structural detail are
        essential for obtaining reliable and biologically meaningful results in cryo-
        electron tomography studies.
    """

    _label = 'resample tomogram'
    _possibleOutputs = {OUTPUT_TOMO_NAME: SetOfTomograms}
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        #     Example 'resample' input:
        # 
        #     resample
        #  
        #  
        #         **   Welcome to Resample   **
        #  
        #              Version : 1.00
        #             Compiled : Dec  2 2017
        #                 Mode : Interactive
        #  
        # Input image file name [tomo2_L1G1_72_d3.rec]       : tomo2_L1G1-dose_filt.rec
        # Output image file name [tomo2_L1G1_72_d3-bin8.rec] : tomo2_L1G1_72_d3-bin8.rec
        # Is the input a volume [YES]                        : YES
        # New X-Size [464]                                   : 464
        # New Y-Size [464]                                   : 464
        # New Z-Size [232]                                   : 232

        # You need a params to belong to a section:
        form.addSection(label=pwutils.Message.LABEL_INPUT)
        form.addParam('inTomograms', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Input tomograms')

        form.addParam('newXsize', IntParam,
                      default=512, validators=[Positive],
                      label='New X-Size',
                      help='Volume will be rescaled to this size in X dimension (voxels)')
        form.addParam('newYsize', IntParam,
                      default=512, validators=[Positive],
                      label='New Y-Size',
                      help='Volume will be rescaled to this size in Y dimension (voxels)')
        form.addParam('newZsize', IntParam,
                      default=256, validators=[Positive],
                      label='New Z-Size',
                      help='Volume will be rescaled to this size in Z dimension (voxels)')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        for tomo in self.inTomograms.get():

            self._insertFunctionStep(self.runTomoResample,
                                     tomo.getFileName(),
                                     needsGPU=False)

        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------

    def runTomoResample(self, tomoFile: str):
        prog = Plugin.getProgram('resample')

        paramDict = {
            'tomoFile': tomoFile,
            'tomoOutName': self._getOutputFn(tomoFile),
            'newXsize': self.newXsize,
            'newYsize': self.newYsize,
            'newZsize': self.newZsize
        }

        # Arguments to the resample command defined in the plugin initialization:
        args = """   << eof
%(tomoFile)s
%(tomoOutName)s
YES
%(newXsize)d
%(newYsize)d
%(newZsize)d
eof\n
"""
        self.runJob(prog, args % paramDict)

    def createOutputStep(self):
        inTomoSet = self.inTomograms.get()
        outTomoSet = self._createSetOfTomograms("resampled")
        outTomoSet.copyInfo(inTomoSet)
        outTomoSet.setSamplingRate(self._getOutputSampling())
        outTomoSet.copyItems(inTomoSet, doClone=False,
                             updateItemCallback=self._updateItem)

        self._defineOutputs(**{OUTPUT_TOMO_NAME: outTomoSet})
        self._defineTransformRelation(self.inTomograms, outTomoSet)

    # --------------------------- INFO functions -----------------------------------
    def _citations(self):
        return ['Grant2018']

    # --------------------------- UTILS functions -------------------------------
    def _updateItem(self, item, row):
        outputFn = self._getOutputFn(item.getFileName())
        if os.path.exists(outputFn):
            item.setFileName(outputFn)
        else:
            item._appendItem = False

    def _getOutputFn(self, tomoFile):
        tomoBaseName = pwutils.removeBaseExt(tomoFile)
        tomoExt = pwutils.getExt(tomoFile)
        output = self._getExtraPath(tomoBaseName + '_resampled' + tomoExt)

        return output

    def _getOutputSampling(self):
        tomos = self.inTomograms.get()
        oldSamplingRate = tomos.getSamplingRate()
        oldXsize = tomos.getFirstItem().getXDim()

        return oldSamplingRate * oldXsize / self.newXsize.get()
