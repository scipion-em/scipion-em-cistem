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
import os

import pyworkflow.utils as pwutils
from pwem.protocols import EMProtocol
from tomo.objects import SetOfTiltSeries
from tomo.protocols import ProtTomoBase

from cistem import Plugin

OUTPUT_TS_NAME = 'resampledTiltSeries'



# The steps return one file here as their result file. On Continue, pyworkflow keeps a finished step
# finished only while its result files exist, so removing the ones of a tilt series that failed makes its
# step run again.
STEP_RESULTS_DIR = 'STEP_RESULTS'

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
        tsFiles = [ts.getFirstItem().getFileName() for ts in self.inputSetOfTiltSeries.get()]
        self._prepareFailedForRetry(tsFiles)

        for tsFile in tsFiles:
            self._insertFunctionStep(self.runTsResample,
                                     tsFile,
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
        try:
            self.runJob(prog, args % paramDict)
        except Exception as e:
            # Not raised: the step would fail and take the whole run with it. The tilt series is left out
            # of the output, and Continue resamples it again.
            self.error("ERROR: Resampling %s has failed: %s" % (tsFile, e))
            pwutils.cleanPath(paramDict['tsOutName'])
        return self._stepResultFile(tsFile)

    def createOutputStep(self):
        inputTs = self.inputSetOfTiltSeries.get()
        outTsSet = self._createSetOfTiltSeries()
        outTsSet.copyInfo(inputTs)
        outTsSet.setSamplingRate(self._getOutputSampling())
        # Only the tilt series that were resampled: a failed one is left out
        outTsSet.copyItems(inputTs, updateTiCallback=self._updateTi,
                           itemSelectedCallback=self._isResampled)

        if not len(outTsSet):
            raise Exception("Resampling failed for every tilt series; see the errors above. "
                            "Continue resamples them again.")

        self._defineOutputs(**{OUTPUT_TS_NAME: outTsSet})
        self._defineTransformRelation(self.inputSetOfTiltSeries, outTsSet)
        return self._stepResultFile('output')

    def _isResampled(self, ts):
        return ts.isEnabled() and os.path.exists(self._getOutputFn(ts.getFirstItem().getFileName()))

    # --------------------------- Continue -----------------------------------
    def _stepResultPath(self, name):
        key = pwutils.removeBaseExt(name) if os.sep in name else name
        return self._getExtraPath(STEP_RESULTS_DIR, '%s.TXT' % key)

    def _stepResultFile(self, name):
        """ The result file of a step: Continue runs the step again once it is removed. """
        pwutils.makePath(self._getExtraPath(STEP_RESULTS_DIR))
        path = self._stepResultPath(name)
        open(path, 'a').close()
        return path

    def _prepareFailedForRetry(self, tsFiles):
        """ Have Continue resample again the tilt series without a result. pyworkflow does not run a finished
        step again, so their steps, and the output one, lose their result files. """
        retried = [tsFile for tsFile in tsFiles if not os.path.exists(self._getOutputFn(tsFile))
                   and os.path.exists(self._stepResultPath(tsFile))]
        for tsFile in retried:
            pwutils.cleanPath(self._stepResultPath(tsFile))
        if retried:
            pwutils.cleanPath(self._stepResultPath('output'))

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
