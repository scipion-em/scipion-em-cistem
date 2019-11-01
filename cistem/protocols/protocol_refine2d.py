# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *  Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca) [2]
# *
# * [1] MRC Laboratory of Molecular Biology (MRC-LMB)
# * [2] Department of Anatomy and Cell Biology, McGill University
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
import re

import pyworkflow.em as em
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        EnumParam, StringParam,
                                        BooleanParam, LabelParam)
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils.path import makePath, cleanPath, createLink

from cistem.convert import (writeReferences)
from cistem.constants import *


class CistemProtRefine2D(em.ProtClassify2D):
    """ Protocol to run 2D classification in cisTEM. """
    _label = 'classify 2D'

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """
        myDict = {
            'run_stack': 'Refine2D/ParticleStacks/particle_stack_%(run)02d.mrc',
            'initial_cls': 'Refine2D/ClassAverages/reference_averages.mrc',
            'iter_cls': 'Refine2D/ClassAverages/class_averages_%(iter)04d.mrc',
            'iter_par': 'Refine2D/Parameters/classification_input_par_%(iter)02d.par',
            'iter_cls_block': 'Refine2D/class_dump_file_%(iter)02d_%(block)02d.dump'
        }
        self._updateFilenamesDict(myDict)

    def _createIterTemplates(self, currRun):
        """ Setup the regex on how to find iterations. """
        clsFn = self._getExtraPath(self._getFileName('classes', run=currRun, iter=1))
        self._iterTemplate = clsFn.replace('classes_01', 'classes_??')
        # Iterations will be identify by classes_XX_ where XX is the iteration
        #  number and is restricted to only 2 digits.
        self._iterRegex = re.compile('classes_(\d{2})')

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('doContinue', BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous '
                           'run of type *%s* class. The refinement will resume '
                           'after the last completed iteration. It is ok to alter '
                           'other parameters.' % self.getClassName())
        form.addParam('continueRun', PointerParam,
                      pointerClass=self.getClassName(),
                      condition='doContinue', allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        form.addParam('inputParticles', PointerParam,
                      label="Input particles",
                      condition='not doContinue',
                      important=True, pointerClass='SetOfParticles',
                      help='Select the input particles.')
        form.addParam('inputClassAvg', PointerParam,
                      condition='not doContinue',
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      label="Input class averages",
                      pointerClass='SetOfAverages',
                      help='Select starting class averages. If not provided, '
                           'they will be generated automatically.')
        form.addParam('numberOfClassAvg', IntParam, default=5,
                      label='Number of classes',
                      help='The number of classes that should be generated. '
                           'This input is only available when starting a '
                           'fresh classification run.')
        form.addParam('numberOfIterations', IntParam, default=20,
                      label='Number of cycles ro run',
                      help='The number of refinement cycles to run. If '
                           'the option "Auto Percent Used" is selected, '
                           '20 cycles are usually sufficient to '
                           'generate good class averages. If '
                           'the user decides to set parameters '
                           'manually, 5 to 10 cycles are usually '
                           'sufficient for a particular set of parameters. '
                           'Several of these shorter runs should be used '
                           'to obtain final class averages, updating '
                           'parameters as needed (e.g. Percent Used, see '
                           'example above).')

        form.addSection(label='Expert options')
        form.addParam('lowResLimit', FloatParam, default=300.0,
                      label='Low resolution limit (A)',
                      help='The data used for classification is usually '
                           'bandpass-limited to exclude spurious '
                           'low-resolution features in the particle '
                           'background. It is therefore good practice '
                           'to set the low-resolution limit to 2.5x '
                           'the approximate particle mask radius.')

        line = form.addLine('High resolution limit (A):',
                            help='The high-resolution bandpass limit should be '
                                 'selected to remove data with low signal-to-noise '
                                 'ratio, and to help speed up the calculation. '
                                 'Since the class averages are not well defined '
                                 'initially, the starting limit should be set to '
                                 'a low resolution, for example 40 A. The limit '
                                 'used for the final iterations should be set '
                                 'sufficiently high, for example 8 A, to include '
                                 'signal originating from protein secondary '
                                 'structure that often helps generate recognizable '
                                 'features in the class averages (see example above).')
        line.addParam('highResLimit1', FloatParam, default=8.0,
                      label='start')
        line.addParam('highResLimit2', FloatParam, default=8.0,
                      label='finish')

        form.addParam('maskRad', FloatParam, default=90.0,
                      label='Mask radius (A)',
                      help='The radius of the circular mask applied to the '
                           'input class averages before classification starts. '
                           'This mask should be sufficiently large to include '
                           'the largest dimension of the particle. The mask '
                           'helps remove noise outside the area of the '
                           'particle.')
        form.addParam('angStep', FloatParam, default=8.0,
                      label='Angular search step (deg)',
                      help='The angular step used to generate the search grid '
                           'when marginalizing over the in-plane rotational '
                           'alignment parameter. The smaller the value, the '
                           'finer the search grid and the slower the search. '
                           'It is often sufficient to set the step to 15deg as '
                           'the algorithm varies the starting point of the '
                           'grid in each refinement cycle, thereby covering '
                           'intermediate in-plane alignment angles. However, '
                           'users can try to reduce the step to 5deg (smaller '
                           'is probably not helpful) to see if class '
                           'averages can be improved further once no further '
                           'improvement is seen at 15deg.')

        line = form.addLine('Search range (A): ',
                            help='The search can be limited in the X and Y '
                                 'directions (measured from the box center) to '
                                 'ensure that only particles close to the box '
                                 'center are used for classification. A '
                                 'smaller range, for example 20 to 40 A, can '
                                 'speed up computation. However, the range '
                                 'should be chosen sufficiently generously to '
                                 'capture most particles. If the range of '
                                 'particle displacements from the box center '
                                 'is unknown, start with a larger value, e.g. '
                                 '100 A, check the results when the run '
                                 'finishes and reduce the range appropriately.')
        line.addParam('rangeX', FloatParam, default=75.0,
                      label='X')
        line.addParam('rangeY', FloatParam, default=75.0,
                      label='Y')

        form.addParam('smooth', FloatParam, default=1.0,
                      label='Smoothing factor [0-1]',
                      help='A factor that reduces the range of likelihoods '
                           'used during classification. A reduced range can '
                           'help prevent the appearance of "empty" classes '
                           '(no members) early in the classification. '
                           'Smoothing may also suppress some high-resolution '
                           'noise. The user should try values between 0.1 '
                           'and 1 if classification suffers from the '
                           'disappearance of small classes or noisy class '
                           'averages.')
        form.addParam('exclExdges', BooleanParam, default=False,
                      label='Exclude blank edges?',
                      help='Should particle boxes with blank edges be '
                           'excluded from classification? Blank edges can '
                           'be the result of particles selected close to '
                           'the edges of micrographs. Blank edges can lead '
                           'to errors in the calculation of the likelihood '
                           'function, which depends on the noise statistics.')
        form.addParam('autoPerc', BooleanParam, default=True,
                      label='Auto percent used?',
                      help='Should the percent of included particles be '
                           'adjusted automatically? A classification '
                           'scheme using initially 300 particles/class, '
                           'then 30% and then 100% is often sufficient to '
                           'obtain good classes and this scheme will be '
                           'used when this option is selected.')
        form.addParam('perc', FloatParam, default=100.0,
                      condition='not autoPerc',
                      label='Percent used',
                      help='The fraction of the dataset used for '
                           'classification. Especially in the beginning, '
                           'classification proceeds more rapidly when only '
                           'a small number of particles are used per class, '
                           'e.g. 300 (see example above). Later runs that '
                           'refine the class averages should use a higher '
                           'percentage and the final run(s) should use all '
                           'the data. This option is only available when '
                           '"Auto Percent Used" is not selected.')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.numberOfBlocks = self._defNumberOfCPUs()
        #self._createFilenameTemplates()
        self._insertContinueStep()
        self._insertItersSteps()
        self._insertFunctionStep("createOutputStep")

    def _insertContinueStep(self):
        if self.doContinue:
            continueRun = self.continueRun.get()
            self.inputParticles.set(None)
            self.numberOfClassAvg.set(continueRun.numberOfClassAvg.get())

            if self.continueIter.get() == 'last':
                self.initIter = continueRun._getCurrIter()
            else:
                self.initIter = int(self.continueIter.get()) + 1
            self._insertFunctionStep('continueStep', self.initIter)

        else:
            self.initIter = 1

        self.finalIter = self.initIter + self.numberOfIterations.get()

    def _insertItersSteps(self):
        """ Insert the steps for all iterations. """
        for iterN in self._allItersN():
            initId = self._insertFunctionStep('initIterStep', iterN)
            paramsDic = self._getParamsIteration(iterN)
            depsRefine = self._insertRefineIterStep(iterN, paramsDic, [initId])
            self._insertFunctionStep("mergeStep", iterN, paramsDic,
                                     prerequisites=depsRefine)

    def _insertRefineIterStep(self, iterN, paramsDic, depsInitId):
        """ Execute the refinement for the current iteration """
        depsRefine = []
        if iterN == 1:
            initAngStepId = self._insertFunctionStep("writeParStep",
                                                     paramsDic,
                                                     prerequisites=depsInitId)
            for block in self._allBlocks():
                refineId = self._insertFunctionStep("refineParticlesStep",
                                                    iterN, block, paramsDic,
                                                    prerequisites=[initAngStepId])
                depsRefine.append(refineId)
        else:
            for block in self._allBlocks():
                refineId = self._insertFunctionStep("refineParticlesStep",
                                                    iterN, block, paramsDic,
                                                    prerequisites=depsInitId)
                depsRefine.append(refineId)
        return depsRefine

    # --------------------------- STEPS functions -----------------------------
    def continueStep(self, iterN):
        """Create a symbolic link of a previous iteration from a previous run."""
        iterN = iterN - 1
        self._setLastIter(iterN)
        continueRun = self.continueRun.get()
        prevDir = continueRun._iterWorkingDir(iterN)
        currDir = self._iterWorkingDir(iterN)
        createLink(prevDir, currDir)

        imgFn = self._getFileName('particles')
        self.writeParticlesByMic(imgFn)

    def initIterStep(self, iterN):
        self._createWorkingDirs()
        inputStack = self._getFileName('run_stack', run=0)
        inputClasses = self._getFileName('initial_cls')

        if iterN == 1:
            imgFn = self._getExtraPath(inputStack)
            self.writeParticlesByMic(imgFn)
            inputCls = self.inputClassAvg.get() if self.inputClassAvg else None

            if inputCls is not None:
                writeReferences(self.inputClassAvg.get(),
                                self._getExtraPath(inputClasses))
        else:
            pass
            #self._splitParFile(iterN, self.numberOfBlocks)

    def writeParStep(self, paramsDic):
        """ Construct a parameter file (.par). """
        #  This function will be called only for iteration 1.
        iterN = 1
        iterDir = self._iterWorkingDir(iterN)
        self._enterDir(iterDir)

        inputParticles = self._getInputParticles()
        magnification = inputParticles.getAcquisition().getMagnification()
        params = {}

        for block in self._allBlocks():
            more = 1
            initPart, lastPart = self._initFinalBlockParticles(block)
            params['initParticle'] = initPart
            params['finalParticle'] = lastPart
            paramDic = self._setParamsRefineParticles(iterN, block)
            paramsRefine = dict(paramsDic.items() + params.items() + paramDic.items())
            f = self.__openParamFile(block, paramsRefine)

            # ToDo: Implement a better method to get the info particles.
            #  Now, you iterate several times over the SetOfParticles
            # (as many threads as you have)
            micIdMap = self._getMicCounter()
            for i, img in self.iterParticlesByMic():
                if img.hasMicId():
                    micId = img.getMicId()
                elif img.hasCoordinate():
                    micId = img.getCoordinate().getMicId()
                else:
                    micId = 0

                film = micIdMap[micId]
                ctf = img.getCTF()
                defocusU, defocusV, astig = ctf.getDefocusU(), ctf.getDefocusV(), ctf.getDefocusAngle()
                partCounter = i + 1

                if partCounter == lastPart:  # The last particle in the block
                    more = 0
                particleLine = ('1, %05d, %05d, %05f, %05f, %02f, %01d\n' %
                                (magnification, film, defocusU, defocusV, astig, more))
                self.__writeParamParticle(f, particleLine)

                if more == 0:  # close the block.
                    self.__closeParamFile(f, paramsRefine)
                    break
        self._leaveDir()

    def refineBlockStep(self, block):
        """ This function execute the bash script for refine a subset(block) of images.
        It will enter in the iteration dir and execute the script there.
        """
        iterDir = self._iterWorkingDir(1)
        program = "./block%03d.sh" % block
        os.chmod(join(iterDir, program), 0775)
        self.runJob(program, "", cwd=iterDir)

    def writeInitialAnglesStep(self):
        """This function write a .par file with all necessary information for a refinement"""

        self.micIdMap = {}
        counter = 0;
        for mic in self._micList:
            self.micIdMap[mic['_micId']] = counter
            counter = counter + 1

        for block in self._allBlocks():
            more = 1
            _, lastPart = self._initFinalBlockParticles(block)
            parFn = self._getFileName('input_par_block', block=block, iter=1, prevIter=0)
            f = open(parFn, 'w')
            f.write("C           PSI   THETA     PHI       SHX       SHY     MAG  FILM      DF1"
                    "      DF2  ANGAST     OCC     -LogP      SIGMA   SCORE  CHANGE\n")

            # ToDo: Implement a better method to get the info particles.
            # Now, you iterate several times over the SetOfParticles (as many threads as you have)
            for i, img in self.iterParticlesByMic():
                partCounter = i + 1

                if partCounter == lastPart:  # The last particle in the block
                    more = 0

                self.writeAnglesLines(partCounter, img, f)

                if more == 0:  # close the block.
                    f.close()
                    break

    def refineParticlesStep(self, iterN, block, paramsDic):
        """Only refine the parameters of the SetOfParticles
        """
        param = {}

        iterDir = self._iterWorkingDir(iterN)
        iniPart, lastPart = self._initFinalBlockParticles(block)
        prevIter = iterN - 1
        param['inputParFn'] = self._getBaseName('input_par_block', block=block, iter=iterN, prevIter=prevIter)
        param['initParticle'] = iniPart
        param['finalParticle'] = lastPart

        paramDic = self._setParamsRefineParticles(iterN, block)

        paramsRefine = dict(paramsDic.items() + paramDic.items() + param.items())
        args = self._prepareCommand()

        if self.mode.get() != 0:
            # frealign program is already in the args script, that's why runJob('')
            self.runJob('', args % paramsRefine, cwd=iterDir)
        else:
            pass
            ##ugly hack when for reconstruction only, just copy the input files
            # inFile  = self._getFileName('input_par_block', block= block, iter=iterN, prevIter=prevIter)
            # outFile = self._getFileName('output_par_block', block=block, iter=iterN)
            # print "I am in dir: ", os.getcwd()
            # print "copying params files", inFile, outFile
            # copyFile(inFile, outFile)

    def mergeStep(self, iterN, paramsDic):
        """Reconstruct a volume from a SetOfParticles with its current parameters refined
        """
        self._mergeAllParFiles(iterN,
                               self.numberOfBlocks)  # merge all parameter files generated in a refineIterStep function.

        initParticle = 1
        finalParticle = self._getInputParticles().getSize()

        os.environ['NCPUS'] = str(self.numberOfBlocks)
        paramsDic['frealign'] = self._getProgram()
        paramsDic['outputParFn'] = self._getBaseName('output_vol_par', iter=iterN)
        paramsDic['initParticle'] = initParticle
        paramsDic['finalParticle'] = finalParticle
        #         paramsDic['paramRefine'] = '0, 0, 0, 0, 0'

        params2 = self._setParams3DR(iterN)

        params3DR = dict(paramsDic.items() + params2.items())

        args = self._prepareCommand()
        iterDir = self._iterWorkingDir(iterN)
        # frealign program is already in the args script, that's why runJob('')
        self.runJob('', args % params3DR, cwd=iterDir)
        self._setLastIter(iterN)

    def createOutputStep(self):
        pass  # should be implemented in subclasses

    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClasses'):
            summary.append("Output classes not ready yet.")
        else:
            summary.append("Input Particles: %s" % self.getObjectTag('inputParticles'))
            summary.append("Classified into *%d* classes." % self.numberOfClassAvg)
            summary.append("Output set: %s" % self.getObjectTag('outputClasses'))

        return summary

    def _methods(self):
        methods = "We classified input particles %s (%d items) " % (
            self.getObjectTag('inputParticles'),
            self._getInputParticles().getSize())
        methods += "into %d classes using refine2d" % self.numberOfClassAvg
        return [methods]

    #--------------------------- UTILS functions ------------------------------
    def _createWorkingDirs(self):
        for dirFn in ['Refine2D/ParticleStacks',
                      'Refine2D/ClassAverages',
                      'Refine2D/Parameters']:
            makePath(self._getExtraPath(dirFn))

    def _allItersN(self):
        """ Iterate over iterations steps """
        for iterN in range(self.initIter, self.finalIter):
            yield iterN

    def _getInputParticlesPointer(self):
        if self.doContinue:
            return self.continueRun.get()._getInputParticlesPointer()
        else:
            return self.inputParticles

    def _getInputParticles(self):
        return self._getInputParticlesPointer().get()

    def _defNumberOfCPUs(self):
        self._micList = []
        self._getMicIdList()
        cpus = max(self.numberOfMpi.get() - 1, self.numberOfThreads.get() - 1, 1)
        numberOfMics = len(self._micList)
        return min(cpus, numberOfMics)

    def _allBlocks(self):
        """ Iterate over all numberOfCPUs. """
        for i in range(1, self.numberOfBlocks + 1):
            yield i

    def writeParticlesByMic(self, stackFn):
        self._getInputParticles().writeStack(stackFn,
                                             orderBy=['_micId', 'id'],
                                             direction='ASC')

    def _getParamsIteration(self, iterN):
        """ Defining the current iteration """
        imgSet = self._getInputParticles()

        # Prepare arguments to call program fralign_v9.exe
        paramsDic = {'mode': self.mode.get(),
                     'innerRadius': self.innerRadius.get(),
                     'outerRadius': self.outerRadius.get(),
                     'molMass': self.molMass.get(),
                     'ThresholdMask': self.ThresholdMask.get(),
                     'pseudoBFactor': self.pseudoBFactor.get(),
                     'avePhaseResidual': self.avePhaseResidual.get(),
                     'angStepSize': self.angStepSize.get(),
                     'numberRandomSearch': int(self.numberRandomSearch.get()),
                     'numberPotentialMatches': int(self.numberPotentialMatches.get()),
                     'sym': self.symmetry.get(),
                     'relMagnification': self.relMagnification.get(),
                     'targetScore': self.targetScore.get(),
                     'score': self.score.get(),
                     'beamTiltX': self.beamTiltX.get(),
                     'beamTiltY': self.beamTiltY.get(),
                     'resol': self.resolution.get(),
                     'lowRes': self.lowResolRefine.get(),
                     'highRes': self.highResolRefine.get(),
                     'resolClass': self.resolClass.get(),
                     'defocusUncertainty': self.defocusUncertainty.get(),
                     'Bfactor': self.Bfactor.get(),
                     'sampling3DR': imgSet.getSamplingRate()
                     }

        # Get the particles stack
        iterDir = self._iterWorkingDir(iterN)
        paramsDic['imageFn'] = os.path.relpath(self._getFileName("particles"), iterDir)
        acquisition = imgSet.getAcquisition()

        # Get the amplitude Contrast of the micrographs
        paramsDic['ampContrast'] = acquisition.getAmplitudeContrast()
        # Get the scanned pixel size of the micrographs
        paramsDic['scannedPixelSize'] = acquisition.getMagnification() * imgSet.getSamplingRate() / 10000
        # Get the voltage and spherical aberration of the microscope
        paramsDic['voltage'] = acquisition.getVoltage()
        paramsDic['sphericalAberration'] = acquisition.getSphericalAberration()

        # Defining the operation mode
        if self.mode == MOD_RECONSTRUCTION:
            paramsDic['mode'] = 0
        elif self.mode == MOD_REFINEMENT:
            paramsDic['mode'] = 1
        elif self.mode == MOD_RANDOM_SEARCH_REFINEMENT:
            paramsDic['mode'] = 2
        elif self.mode == MOD_SIMPLE_SEARCH_REFINEMENT:
            paramsDic['mode'] = 3
        else:
            paramsDic['mode'] = 4

        # Defining the operation mode for iteration 1.
        if self.Firstmode == MOD2_SIMPLE_SEARCH_REFINEMENT:
            paramsDic['mode2'] = -3
        else:
            paramsDic['mode2'] = -4

        # Defining if magnification refinement is going to do
        if self.doMagRefinement and iterN != 1:
            paramsDic['doMagRefinement'] = 'T'
        else:
            paramsDic['doMagRefinement'] = 'F'

        # Defining if defocus refinement is going to do
        if self.doDefRefinement and iterN != 1:
            paramsDic['doDefocusRef'] = 'T'
        else:
            paramsDic['doDefocusRef'] = 'F'

        # Defining if astigmatism refinement is going to do
        if self.doAstigRefinement and iterN != 1:
            paramsDic['doAstigRef'] = 'T'
        else:
            paramsDic['doAstigRef'] = 'F'

        # Defining if defocus refinement for individual particles is going to do
        if self.doDefPartRefinement and iterN != 1:
            paramsDic['doDefocusPartRef'] = 'T'
        else:
            paramsDic['doDefocusPartRef'] = 'F'

        if self.methodEwaldSphere == EWA_DISABLE:
            paramsDic['metEwaldSphere'] = 0
        elif self.methodEwaldSphere == EWA_SIMPLE:
            paramsDic['metEwaldSphere'] = 1
        elif self.methodEwaldSphere == EWA_REFERENCE:
            paramsDic['metEwaldSphere'] = 2
        elif self.methodEwaldSphere == EWA_SIMPLE_HAND:
            paramsDic['metEwaldSphere'] = -1
        else:
            paramsDic['metEwaldSphere'] = -2

        # Defining if apply extra real space symmetry
        if self.doExtraRealSpaceSym:
            paramsDic['doExtraRealSpaceSym'] = 'T'
        else:
            paramsDic['doExtraRealSpaceSym'] = 'F'

        # Defining if wiener filter is going to apply
        if self.doWienerFilter:
            paramsDic['doWienerFilter'] = 'T'
        else:
            paramsDic['doWienerFilter'] = 'F'

        # Defining if wiener filter is going to calculate and apply
        if self.doBfactor:
            paramsDic['doBfactor'] = 'T'
        else:
            paramsDic['doBfactor'] = 'F'

        # Defining if matching projections is going to write
        if self.writeMatchProjections:
            paramsDic['writeMatchProj'] = 'T'
        else:
            paramsDic['writeMatchProj'] = 'F'

        # Defining the method to FSC calcutalion
        if self.methodCalcFsc == FSC_CALC:
            paramsDic['metFsc'] = 0
        elif self.methodCalcFsc == FSC_3DR_ODD:
            paramsDic['metFsc'] = 1
        elif self.methodCalcFsc == FSC_3DR_EVEN:
            paramsDic['metFsc'] = 2
        elif self.methodCalcFsc == FSC_3DR_ALL:
            paramsDic['metFsc'] = 3

        if self.doAditionalStatisFSC:
            paramsDic['doAditionalStatisFSC'] = 'T'
        else:
            paramsDic['doAditionalStatisFSC'] = 'F'

        if self.memory == MEM_0:
            paramsDic['memory'] = 0
        elif self.memory == MEM_1:
            paramsDic['memory'] = 1
        elif self.memory == MEM_2:
            paramsDic['memory'] = 2
        else:
            paramsDic['memory'] = 3

        if self.interpolationScheme == INTERPOLATION_0:
            paramsDic['interpolation'] = 0
        else:
            paramsDic['interpolation'] = 1

        if self.paramRefine == REF_ALL:
            paramsDic['paramRefine'] = '1, 1, 1, 1, 1'
        elif self.paramRefine == REF_ANGLES:
            paramsDic['paramRefine'] = '1, 1, 1, 0, 0'
        elif self.paramRefine == REF_SHIFTS:
            paramsDic['paramRefine'] = '0, 0, 0, 1, 1'
        else:
            paramsDic['paramRefine'] = '0, 0, 0, 0, 0'

        return paramsDic
