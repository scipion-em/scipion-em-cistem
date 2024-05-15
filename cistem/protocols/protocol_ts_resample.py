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
    """
    Resample tilt series by Fourier cropping/padding using cisTEM.
    This is equivalent to binning/unbinning operations but free of aliasing artifacts.

    More info:
        https://cistem.org
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
                                     ts.getFirstItem().getFileName())

        self._insertFunctionStep(self.createOutputStep)

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
