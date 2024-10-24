# **************************************************************************
# *
# * Authors:    Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *             Scipion Team (scipion@cnb.csic.es) [2]
# *
# * [1] MRC Laboratory of Molecular Biology (MRC-LMB)
# * [2] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from pyworkflow.utils import weakImport
with weakImport('tomo'):
    from pyworkflow.tests import DataSet, setupTestProject
    from pyworkflow.utils import magentaStr

    from tomo.protocols import ProtImportTs, ProtImportTsCTF
    from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
    from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer
    from ..protocols import CistemProtTsCtffind


    class TestBase(TestBaseCentralizedLayer):
        @classmethod
        def setUpClass(cls):
            setupTestProject(cls)
            cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
            cls.inputSoTS = cls._runImportTs()

        @classmethod
        def _runImportTs(cls):
            print(magentaStr("\n==> Importing data - tilt-series:"))
            protImportTs = cls.newProtocol(ProtImportTs,
                                           filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                           filesPattern=DataSetRe4STATuto.tsPattern.value,
                                           exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value,
                                           anglesFrom=2,  # From tlt file
                                           voltage=DataSetRe4STATuto.voltage.value,
                                           sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                           amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                           samplingRate=DataSetRe4STATuto.unbinnedPixSize.value,
                                           doseInitial=DataSetRe4STATuto.initialDose.value,
                                           dosePerFrame=DataSetRe4STATuto.dosePerTiltImg.value,
                                           tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

            cls.launchProtocol(protImportTs)
            tsImported = getattr(protImportTs, 'outputTiltSeries', None)
            return tsImported


    class TestCtffindTs(TestBase):
        def testCistemCtfFind(self):
            print(magentaStr("\n==> Testing cistem - ctffind:"))
            protCTF = CistemProtTsCtffind(inputTiltSeries=self.inputSoTS,
                                          lowRes=30,
                                          highRes=3,
                                          minDefocus=10000,
                                          maxDefocus=20000,
                                          numberOfThreads=3)  # 1 per each TS plus the main one
            self.launchProtocol(protCTF, wait=True)
            outCtfs = getattr(protCTF, protCTF._possibleOutputs.CTFs.name, None)
            self.assertIsNotNone(outCtfs, "SetOfCTFTomoSeries has not been produced.")
            self.checkCTFs(outCtfs, expectedSetSize=2)

        def testCistemImportCtfFiles(self):
            print(magentaStr("\n==> Importing data - ctffind files:"))
            protImport = ProtImportTsCTF(filesPath=self.ds.getFile(DataSetRe4STATuto.cistemFilesPath.name),
                                         filesPattern='*.txt',
                                         importFrom=0,  # ctffind
                                         inputSetOfTiltSeries=self.inputSoTS)
            self.launchProtocol(protImport)
            outCtfs = getattr(protImport, protImport._possibleOutputs.CTFs.name, None)
            self.assertIsNotNone(outCtfs, "SetOfCTFTomoSeries has not been produced.")
            self.checkCTFs(outCtfs, expectedSetSize=2)
