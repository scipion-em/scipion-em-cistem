# **************************************************************************
# *
# * Authors:     Ricardo D. Righetto (ricardo.righetto@unibas.ch)
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


'''
A protocol to resample tilt series by Fourier cropping/padding using cisTEM.
'''
# Inspired by protocol_resizeTS.py from scipion-em-xmipptomo:
# https://github.com/I2PC/scipion-em-xmipptomo/blob/devel/xmipptomo/protocols/protocol_resizeTS.py

from cistem import Plugin
from tomo.protocols import ProtTomoBase
import tomo.objects as tomoObj
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, IntParam
from pyworkflow.utils import *
from pyworkflow.object import Set
from tomo.objects import SetOfTiltSeries, TiltSeries
from pwem.emlib.image import ImageHandler

OUTPUT_TS_NAME = 'resampledTiltSeries'
OUTPUT_DIR = 'extra/'

class ProtTsResample(EMProtocol, ProtTomoBase):
    '''
    Resample tilt series by Fourier cropping/padding using cisTEM. This is equivalent to binning/unbinning operations but free of aliasing artifacts.

    More info:
        https://cistem.org
    '''

    _label = 'resample tilt series'
    _possibleOutputs = {OUTPUT_TS_NAME: SetOfTiltSeries}
    _devStatus = BETA

    tsList = []

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ''' Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        '''
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
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfTiltSeries', PointerParam,
                        pointerClass='SetOfTiltSeries',
                        allowsNull=False,
                        label='Input tilt series')

        form.addParam('newXsize', IntParam,
                        default=512,
                        label='New X-Size',
                        help='Volume will be rescaled to this size in X dimension (voxels)')
        form.addParam('newYsize', IntParam,
                        default=512,
                        label='New Y-Size',
                        help='Volume will be rescaled to this size in Y dimension (voxels)')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):

        for ts in self.inputSetOfTiltSeries.get():

            self._insertFunctionStep(self.runTsResample,
                        ts.getObjId())

            self._insertFunctionStep(self.createOutputStep,
                        ts.getObjId())
        
        self._insertFunctionStep('closeStreamStep')

    def runTsResample(self, tsObjId):

        prog = Plugin.getProgram('resample')

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        firstItem = ts.getFirstItem()
        tsFile = firstItem.getFileName()

        tsBaseName = removeBaseExt(tsFile)
        tsExt = getExt(tsFile)
        paramDict = {}
        paramDict['tsFile'] = tsFile
        paramDict['tsOutName'] = self.getWorkingDir() + '/' + OUTPUT_DIR + '/' + \
            tsBaseName + '_resampled' + tsExt
        
        paramDict['newXsize'] = self.newXsize
        paramDict['newYsize'] = self.newYsize

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

        self.tsList.append(paramDict['tsOutName'])

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]

        firstItem = ts.getFirstItem()
        tsFile = firstItem.getFileName()
        tsBaseName = removeBaseExt(tsFile)
        tsExt = getExt(tsFile)

        tsOutName = self.getWorkingDir() + '/' + OUTPUT_DIR + '/' + \
            tsBaseName + '_resampled' + tsExt

        tsId = ts.getTsId()

        oldSamplingRate = ts.getSamplingRate()
        oldXsize = ts.getDim()[0]
        self.samplingRate = oldSamplingRate * oldXsize / float(self.newXsize)

        outputSetOfTiltSeries = self.getOutputSetOfTiltSeries()

        newTs = TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputSetOfTiltSeries.append(newTs)

        newTs.setSamplingRate(self.samplingRate)

        for index, ti in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(ti, copyId=True)
            newTi.setLocation(index + 1, tsOutName)

            if ti.hasTransform():
                newTi.setTransform(ti.getTransform())

            newTi.setSamplingRate(self.samplingRate)

            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))
        newTs.write(properties=False)

        outputSetOfTiltSeries.update(newTs)
        outputSetOfTiltSeries.updateDim()
        outputSetOfTiltSeries.write()

        self._store()

    def closeStreamStep(self):
        self.getOutputSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions -------------------------------
    def getOutputSetOfTiltSeries(self):
        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            outputSetOfTiltSeries.setSamplingRate(self.samplingRate)
            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)
        return self.outputSetOfTiltSeries

    # --------------------------- INFO functions -----------------------------------
    def _citations(self):

        cites = ['Grant2018']

        return cites

    def _validate(self):
        errors = []
        if self.newXsize <= 0:
            errors.append('New X size must be greater than zero!')
        if self.newYsize <= 0:
            errors.append('New Y size must be greater than zero!')
        return errors
