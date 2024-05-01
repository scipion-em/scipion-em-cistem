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
    """
    Resample tomograms by Fourier cropping/padding using cisTEM.
    This is equivalent to binning/unbinning operations but free of aliasing artifacts.

    More info:
        https://cistem.org
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
                                     tomo.getFileName())

        self._insertFunctionStep(self.createOutputStep)

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
