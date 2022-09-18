# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [2]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [3]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Centro Nacional de Biotecnologia, CSIC, Spain
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
from enum import Enum

from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pwem.objects import CTFModel, Set

from ..convert import readCtfModelStack, parseCtffind4Output

from tomo.objects import CTFTomo, CTFTomoSeries, SetOfCTFTomoSeries
from tomo.protocols.protocol_base import ProtTomoImportFiles


class outputs(Enum):
    outputCTFs = SetOfCTFTomoSeries


class CistemProtTsImportCtf(ProtTomoImportFiles):
    """ Protocol to import CTF estimation of a tilt-series from CTFFIND4. """
    _label = 'import tomo CTFs'
    _devStatus = BETA
    _possibleOutputs = outputs

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineImportParams(form)

        form.addParam('exclusionWords', params.StringParam,
                      label='Exclusion words:',
                      help="List of words separated by a space that "
                           "the path should not have",
                      expertLevel=params.LEVEL_ADVANCED)

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input tilt-series',
                      help='Select the tilt-series to which the imported '
                           'estimation of CTF will be paired. The file '
                           'names of the file and the defocus file must '
                           'be the same (except the extension).')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.importSetOfCtfTomoSeries)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def importSetOfCtfTomoSeries(self):
        outputCtfs = self._getOutputSet()
        if outputCtfs is None:
            outputCtfs = self._createOutputSet()
        else:
            outputCtfs.enableAppend()

        for ts in self._getInputTs():
            tsId = ts.getTsId()
            tsObjId = ts.getObjId()
            tsFileName = ts.getFirstItem().parseFileName(extension='')

            for defocusFn in self.iterFiles():
                if tsFileName == pwutils.removeBaseExt(defocusFn):
                    print("Parsing file: " + defocusFn)
                    newCTFTomoSeries = CTFTomoSeries()
                    newCTFTomoSeries.copyInfo(ts)
                    newCTFTomoSeries.setTiltSeries(ts)
                    newCTFTomoSeries.setObjId(tsObjId)
                    newCTFTomoSeries.setTsId(tsId)

                    outputCtfs.append(newCTFTomoSeries)

                    outputPsd = self._findPsdFile(defocusFn)
                    ctfResult = parseCtffind4Output(defocusFn)
                    ctf = CTFModel()

                    for i, ti in enumerate(ts):
                        newCtfTomo = self.getCtfTi(ctf, ctfResult, i, outputPsd)
                        newCtfTomo.setIndex(i+1)
                        newCTFTomoSeries.append(newCtfTomo)

                    newCTFTomoSeries.calculateDefocusUDeviation()
                    newCTFTomoSeries.calculateDefocusVDeviation()

                    if not (newCTFTomoSeries.getIsDefocusUDeviationInRange() and
                            newCTFTomoSeries.getIsDefocusVDeviationInRange()):
                        newCTFTomoSeries.setEnabled(False)

                    newCTFTomoSeries.write(properties=False)
                    outputCtfs.update(newCTFTomoSeries)

    def closeOutputSetsStep(self):
        self._getOutputSet().setStreamState(Set.STREAM_CLOSED)
        self._getOutputSet().write()
        self._store()

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self._getOutputSet() is not None:
            summary.append("Imported CTF tomo series: %d.\n"
                           % (self._getOutputSet().getSize()))
        else:
            summary.append("Output CTFs not ready yet.")
        return summary

    def _validate(self):
        errorMsg = []
        if not self.getMatchFiles():
            errorMsg.append('Unable to find the files provided:\n\n'
                            '\t-filePath = %s\n'
                            '\t-pattern = %s\n' % (self.filesPath.get(),
                                                   self.filesPattern.get()))

        return errorMsg

    # --------------------------- UTILS functions -----------------------------
    def _getInputTs(self, pointer=False):
        return (self.inputSetOfTiltSeries.get() if not pointer
                else self.inputSetOfTiltSeries)

    def getCtfTi(self, ctf, ctfArray, tiIndex, psdStack=None):
        """ Parse the CTF object estimated for this Tilt-Image. """
        readCtfModelStack(ctf, ctfArray, item=tiIndex)
        if psdStack is not None:
            ctf.setPsdFile(f"{tiIndex+1}@" + psdStack)
        ctfTomo = CTFTomo.ctfModelToCtfTomo(ctf)

        return ctfTomo

    def _getOutputSet(self):
        return getattr(self, self._possibleOutputs.outputCTFs.name, None)

    def _createOutputSet(self):
        outputCtfs = SetOfCTFTomoSeries.create(self._getPath(),
                                               template='CTFmodels%s.sqlite')
        outputCtfs.setSetOfTiltSeries(self._getInputTs(pointer=True))
        outputCtfs.setStreamState(Set.STREAM_OPEN)
        self._defineOutputs(**{outputs.outputCTFs.name: outputCtfs})
        return outputCtfs

    def iterFiles(self):
        """ Iterate through the files matched with the pattern.
        Provide the fileName and fileId.
        """
        filePaths = self.getMatchFiles()
        filePaths = self._excludeByWords(filePaths)

        for fileName in filePaths:
            if self._idRegex:
                # Try to match the file id from filename
                # this is set by the user by using #### format in the pattern
                match = self._idRegex.match(fileName)
                if match is None:
                    raise Exception("File '%s' doesn't match the pattern '%s'"
                                    % (fileName, self.getPattern()))

            yield fileName

    def _excludeByWords(self, files):
        exclusionWords = self.exclusionWords.get()

        if exclusionWords is None:
            return files

        exclusionWordList = exclusionWords.split()
        allowedFiles = []

        for file in files:
            if any(bannedWord in file for bannedWord in exclusionWordList):
                print("%s excluded. Contains any of %s" % (file, exclusionWords))
                continue
            allowedFiles.append(file)

        return allowedFiles

    def _findPsdFile(self, fn):
        fnBase = pwutils.removeExt(fn)
        for suffix in ['_psd.mrc', '.mrc', '_ctf.mrcs',
                       '.mrcs', '.ctf']:
            psdFile = fnBase + suffix
            if os.path.exists(psdFile):
                if psdFile.endswith('.ctf'):
                    psdFile += ':mrc'
                return psdFile
        return None
