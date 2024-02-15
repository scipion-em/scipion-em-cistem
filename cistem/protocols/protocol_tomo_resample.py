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
A protocol to resample tomograms by Fourier cropping/padding using cisTEM.
'''
from cistem import Plugin
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, IntParam
from pyworkflow.utils import *
from tomo.objects import SetOfTomograms, Tomogram
from pwem.convert.headers import Ccp4Header

OUTPUT_TOMO_NAME = 'resampledTomos'
OUTPUT_DIR = 'extra/'

class ProtTomoResample(EMProtocol):
    '''
    Resample tomograms by Fourier cropping/padding using cisTEM. This is equivalent to binning/unbinning operations but free of aliasing artifacts.

    More info:
        https://cistem.org
    '''

    _label = 'resample tomogram'
    _possibleOutputs = {OUTPUT_TOMO_NAME: SetOfTomograms}
    _devStatus = BETA

    tomoList = []

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
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inTomograms', PointerParam,
                        pointerClass='SetOfTomograms',
                        allowsNull=False,
                        label='Input tomograms')

        form.addParam('newXsize', IntParam,
                        default=512,
                        label='New X-Size',
                        help='Volume will be rescaled to this size in X dimension (voxels)')
        form.addParam('newYsize', IntParam,
                        default=512,
                        label='New Y-Size',
                        help='Volume will be rescaled to this size in Y dimension (voxels)')
        form.addParam('newZsize', IntParam,
                        default=256,
                        label='New Z-Size',
                        help='Volume will be rescaled to this size in Z dimension (voxels)')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):

        for tomo in self.inTomograms.get():

            self._insertFunctionStep(self.runTomoResample,
                        tomo.getFileName())

        self._insertFunctionStep(self.createOutputStep)

    def runTomoResample(self, tomoFile: str):

        prog = Plugin.getProgram('resample')


        tomoBaseName = removeBaseExt(tomoFile)
        tomoExt = getExt(tomoFile)
        paramDict = {}
        paramDict['tomoFile'] = tomoFile
        paramDict['tomoOutName'] = self.getWorkingDir() + '/' + OUTPUT_DIR + '/' + \
            tomoBaseName + '_resampled' + tomoExt
        
        paramDict['newXsize'] = self.newXsize
        paramDict['newYsize'] = self.newYsize
        paramDict['newZsize'] = self.newZsize

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

        self.tomoList.append(paramDict['tomoOutName'])

    def createOutputStep(self):
        labelledSet = self._genOutputSetOfTomograms(
            self.tomoList, 'resampled')
        self._defineOutputs(**{OUTPUT_TOMO_NAME: labelledSet})
        self._defineSourceRelation(self.inTomograms.get(), labelledSet)

    def _genOutputSetOfTomograms(self, tomoList, suffix):
        tomoSet = SetOfTomograms.create(
            self._getPath(), template='tomograms%s.sqlite', suffix=suffix)
        inTomoSet = self.inTomograms.get()
        tomoSet.copyInfo(inTomoSet)

        outSamplingRate = Ccp4Header(tomoList[0], readHeader=True).getSampling()
        tomoSet.setSamplingRate(outSamplingRate[0])

        counter = 1
        for file, inTomo in zip(tomoList, inTomoSet):
            tomo = Tomogram()
            tomo.copyInfo(inTomo)
            tomo.setLocation((counter, file))
            tomoSet.append(tomo)
            counter += 1    

        return tomoSet

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
        if self.newZsize <= 0:
            errors.append('New Z size must be greater than zero!')
        return errors