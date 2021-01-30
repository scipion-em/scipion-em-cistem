# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es) [1]
# * Authors:    Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] MRC Laboratory of Molecular Biology (MRC-LMB)
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

from pyworkflow.tests import BaseTest, DataSet, setupTestProject

from pwem.protocols import ProtImportMicrographs, ProtImportParticles
from pyworkflow.utils import magentaStr

from ..protocols import (CistemProtCTFFind, CistemProtFindParticles,
                         CistemProtRefine2D)


class TestBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.micFn = cls.dataset.getFile('allMics')

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage,
                            sphericalAberration):
        """ Run an Import micrograph protocol. """
        kwargs = {
            'filesPath': pattern,
            'voltage': voltage,
            'sphericalAberration': sphericalAberration,
            'samplingRate': samplingRate
        }

        cls.protImport = ProtImportMicrographs(**kwargs)
        cls.launchProtocol(cls.protImport)
        cls.assertIsNotNone(cls.protImport.outputMicrographs,
                            "SetOfMicrographs has not been produced.")

        return cls.protImport

    @classmethod
    def runImportMicrographBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern,
                                       samplingRate=1.237,
                                       voltage=300,
                                       sphericalAberration=2)

    @classmethod
    def runImportParticlesSqlite(cls, sqliteFn, sampling):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportParticles,
                                     importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                     sqliteFile=sqliteFn,
                                     samplingRate=sampling,
                                     haveDataBeenPhaseFlipped=False)
        cls.launchProtocol(protImport)
        cls.assertIsNotNone(protImport.outputParticles,
                            "SetOfParticles has not been produced.")

        return protImport

    @classmethod
    def runCtffind(cls, inputMics):
        """ Run CTFFind protocol. """
        cls.protCTF = CistemProtCTFFind()
        cls.protCTF.inputMicrographs.set(inputMics)
        cls.launchProtocol(cls.protCTF)
        cls.assertIsNotNone(cls.protCTF.outputCTF,
                            "SetOfCTF has not been produced.")

        return cls.protCTF


class TestCtffind4(TestBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestBase.setData()
        print(magentaStr("\n==> Importing data - micrographs:"))
        cls.protImport = cls.runImportMicrographBPV(cls.micFn)
    
    def testCtffind4V1(self):
        print(magentaStr("\n==> Testing cistem - ctffind:"))
        protCTF = CistemProtCTFFind()
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        self.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "SetOfCTF has not been produced.")

        valuesList = [[24067, 23587], [22373, 22039], [22653, 22480]]
        for ctfModel, values in zip(protCTF.outputCTF, valuesList):
            self.assertAlmostEquals(ctfModel.getDefocusU(), values[0], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusV(), values[1], delta=1000)


class TestFindParticles(TestBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestBase.setData()
        print(magentaStr("\n==> Importing data - micrographs:"))
        cls.protImport = cls.runImportMicrographBPV(cls.micFn)
        print(magentaStr("\n==> Running cistem - ctffind:"))
        cls.protCtfRun = cls.runCtffind(cls.protImport.outputMicrographs)

    def testFindParts(self):
        print(magentaStr("\n==> Testing cistem - find particles:"))
        protPick = self.newProtocol(CistemProtFindParticles,
                                    maxradius=300.,
                                    radius=300.,
                                    threshold=3.0,
                                    ptclWhite=True,  # wrong but it works!
                                    minDist=250)
        protPick.inputMicrographs.set(self.protImport.outputMicrographs)
        protPick.ctfRelations.set(self.protCtfRun.outputCTF)

        self.launchProtocol(protPick, wait=True)
        self.assertIsNotNone(protPick.outputCoordinates,
                             "SetOfCoordinates has not been produced.")


class TestRefine2D(TestBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.particlesFn = cls.dataset.getFile('import/case2/particles.sqlite')
        print(magentaStr("\n==> Importing data - particles:"))
        cls.protImport = cls.runImportParticlesSqlite(cls.particlesFn,
                                                      sampling=3.5)

    def testClassify2D(self):
        print(magentaStr("\n==> Testing cistem - refine2d:"))
        prot2D = self.newProtocol(CistemProtRefine2D,
                                  numberOfClassAvg=4,
                                  numberOfIterations=3,
                                  areParticlesBlack=False,
                                  numberOfThreads=3)
        prot2D.inputParticles.set(self.protImport.outputParticles)

        self.launchProtocol(prot2D)
        self.assertIsNotNone(prot2D.outputClasses,
                             "SetOfClasses2D has not been produced.")
