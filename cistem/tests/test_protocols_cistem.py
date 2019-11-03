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

from pyworkflow.em import *
from pyworkflow.tests import *

from cistem import *
from cistem.protocols import (CistemProtCTFFind, CistemProtFindParticles,
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

        # Check that input micrographs have been imported
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. '
                            'outputMicrographs is None.' % pattern)

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

        if protImport.outputParticles is None:
            raise Exception('Import of particles: %s, failed. '
                            'outputParticles is None.' % sqliteFn)

        return protImport

    @classmethod
    def runCtffind(cls, inputMics):
        """ Run CTFFind protocol. """
        cls.protCTF = CistemProtCTFFind()
        cls.protCTF.inputMicrographs.set(inputMics)
        cls.launchProtocol(cls.protCTF)

        if cls.protCTF.outputCTF is None:
            raise Exception("SetOfCTF has not been produced.")

        return cls.protCTF


class TestCtffind4(TestBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestBase.setData()
        cls.protImport = cls.runImportMicrographBPV(cls.micFn)
    
    def testCtffind4V1(self):
        protCTF = CistemProtCTFFind()
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        self.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "SetOfCTF has not been produced.")

        valuesList = [[24067, 23587, 58], [22373, 22039, 66], [22653, 22480, 5]]
        for ctfModel, values in izip(protCTF.outputCTF, valuesList):
            self.assertAlmostEquals(ctfModel.getDefocusU(),values[0], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusV(),values[1], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusAngle(),values[2], delta=10)


class TestFindParticles(TestBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestBase.setData()
        cls.protImport = cls.runImportMicrographBPV(cls.micFn)
        cls.protCtfRun = cls.runCtffind(cls.protImport.outputMicrographs)

    def testFindParts(self):
        protPick = CistemProtFindParticles()
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
        cls.protImport = cls.runImportParticlesSqlite(cls.particlesFn,
                                                      sampling=3.5)

    def testClassify2D(self):
        prot2D = self.newProtocol(CistemProtRefine2D,
                                  numberOfThreads=3)
        prot2D.numberOfClassAvg.set(4)
        prot2D.numberOfIterations.set(3)
        prot2D.areParticlesBlack.set(False)
        prot2D.inputParticles.set(self.protImport.outputParticles)

        self.launchProtocol(prot2D)
        self.assertIsNotNone(prot2D.outputClasses,
                             "SetOfClasses2D has not been produced.")
