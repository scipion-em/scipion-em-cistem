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
import logging
import os
from enum import Enum
from os.path import dirname, join
from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pwem.objects import CTFModel
from ..convert import readCtfModelStack, parseCtffind4Output
from tomo.objects import CTFTomo, CTFTomoSeries, SetOfCTFTomoSeries, TiltSeries
from tomo.protocols.protocol_base import ProtTomoImportFiles

logger = logging.getLogger(__name__)


class outputs(Enum):
    CTFs = SetOfCTFTomoSeries


class CistemProtTsImportCtf(ProtTomoImportFiles):
    """ Protocol to import CTF estimation of a tilt-series from CTFFIND4. """
    _label = 'import tomo CTFs'
    _devStatus = PROD
    _possibleOutputs = outputs

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)

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

    # --------------------------- STEPS functions -----------------------------
    def importSetOfCtfTomoSeries(self):
        outputCtfs = self._getOutputSet()
        if outputCtfs is None:
            outputCtfs = self._createOutputSet()

        defocusFiles = [fn for fn in self.iterFiles()]
        for ts in self._getInputTs():
            tsId = ts.getTsId()
            tsObjId = ts.getObjId()

            for defocusFn in defocusFiles:
                if tsId == pwutils.removeBaseExt(defocusFn).replace('_ctf', ''):
                    logger.info("Parsing file: " + defocusFn)
                    newCTFTomoSeries = CTFTomoSeries()
                    newCTFTomoSeries.copyInfo(ts)
                    newCTFTomoSeries.setTiltSeries(ts)
                    newCTFTomoSeries.setObjId(tsObjId)
                    newCTFTomoSeries.setTsId(tsId)

                    outputCtfs.append(newCTFTomoSeries)

                    outputPsd = self._findPsdFile(defocusFn, tsId)
                    ctfResult = parseCtffind4Output(defocusFn)
                    ctf = CTFModel()

                    for i, ti in enumerate(ts):
                        newCtfTomo = self.getCtfTi(ctf, ctfResult, i, outputPsd)
                        newCtfTomo.setIndex(i + 1)
                        newCTFTomoSeries.append(newCtfTomo)

                    newCTFTomoSeries.calculateDefocusUDeviation()
                    newCTFTomoSeries.calculateDefocusVDeviation()

                    if not (newCTFTomoSeries.getIsDefocusUDeviationInRange() and
                            newCTFTomoSeries.getIsDefocusVDeviationInRange()):
                        newCTFTomoSeries.setEnabled(False)

                    outputCtfs.update(newCTFTomoSeries)
                    defocusFiles.remove(defocusFn)
                    break

        outputCtfs.write()
        self._store()

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        outSet = self._getOutputSet()
        if outSet:
            summary.append("Imported CTF tomo series: *%i*" % (self._getOutputSet().getSize()))
            inTsIdSet = set(self._getInputTs().getUniqueValues(TiltSeries.TS_ID_FIELD))
            outTsIdSet = set(outSet.getUniqueValues(TiltSeries.TS_ID_FIELD))
            nonCommonOutTsIds = inTsIdSet.difference(outTsIdSet)
            if nonCommonOutTsIds:
                summary.append(f'Some input TS did not match to the introduced files:\n'
                               f'\t*Non-matching tsIds = {nonCommonOutTsIds}*')
        else:
            summary.append("Output CTFs not ready yet.")
        return summary

    def _validate(self):
        errorMsg = []
        matchingFiles = self.getMatchFiles()
        if matchingFiles:
            tsIdList = self._getInputTs().getUniqueValues(TiltSeries.TS_ID_FIELD)
            defocusBNames = [pwutils.removeBaseExt(defocusFn).replace('_ctf', '')
                             for defocusFn in self.iterFiles()]
            matchResults = list(set(tsIdList) & set(defocusBNames))
            if not matchResults:
                errorMsg.append(f'No matching files found.\n'
                                f'Ctffind files are expected to be named like tsId_ctf_avrot.txt.\n'
                                f'Present tilt series identifiers (tsId) found in the introduced tilt series are '
                                f'{tsIdList}.\n'
                                f'The suffixes "_ctf" or "_avrot" may not be present, they are not mandatory.')
        else:
            errorMsg.append('Unable to find the files provided:\n\n'
                            '\t-filePath = %s\n'
                            '\t-pattern = %s\n' % (self.filesPath.get(),
                                                   self.filesPattern.get()))

        return errorMsg

    def allowsDelete(self, obj):
        return True

    # --------------------------- UTILS functions -----------------------------
    def _getInputTs(self, pointer=False):
        return (self.inputSetOfTiltSeries.get() if not pointer
                else self.inputSetOfTiltSeries)

    @staticmethod
    def getCtfTi(ctf, ctfArray, tiIndex, psdStack=None):
        """ Parse the CTF object estimated for this Tilt-Image. """
        readCtfModelStack(ctf, ctfArray, item=tiIndex)
        if psdStack is not None:
            ctf.setPsdFile(f"{tiIndex + 1}@" + psdStack)
        ctfTomo = CTFTomo.ctfModelToCtfTomo(ctf)

        return ctfTomo

    def _getOutputSet(self):
        return getattr(self, self._possibleOutputs.CTFs.name, None)

    def _createOutputSet(self):
        outputCtfs = SetOfCTFTomoSeries.create(self._getPath(),
                                               template='CTFmodels%s.sqlite')
        outputCtfs.setSetOfTiltSeries(self._getInputTs(pointer=True))
        self._defineOutputs(**{outputs.CTFs.name: outputCtfs})
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
                    raise ValueError("File '%s' doesn't match the pattern '%s'"
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

    @staticmethod
    def _findPsdFile(fn, tsId):
        fnBase = join(dirname(fn), tsId)
        for suffix in ['_psd.mrc', '.mrc', '_ctf.mrcs',
                       '.mrcs', '.ctf']:
            psdFile = fnBase + suffix
            if os.path.exists(psdFile):
                if psdFile.endswith('.ctf'):
                    psdFile += ':mrc'
                return psdFile
        return None
