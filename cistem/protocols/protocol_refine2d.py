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

import re
from glob import glob
from collections import OrderedDict

import pyworkflow.em as em
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        IntParam, BooleanParam,)
from pyworkflow.utils.path import (makePath, createLink,
                                   cleanPattern, moveFile,
                                   exists)

from cistem import Plugin
from cistem.convert import (writeReferences, geometryFromMatrix,
                            rowToAlignment, HEADER_COLUMNS)
from cistem.constants import *


class CistemProtRefine2D(em.ProtClassify2D):
    """ Protocol to run 2D classification in cisTEM. """
    _label = 'classify 2D'

    def __init__(self, **args):
        em.ProtClassify2D.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """
        myDict = {
            'run_stack': 'Refine2D/ParticleStacks/particle_stack_%(run)02d.mrc',
            'initial_cls': 'Refine2D/ClassAverages/reference_averages.mrc',
            'iter_cls': 'Refine2D/ClassAverages/class_averages_%(iter)04d.mrc',
            'iter_par': 'Refine2D/Parameters/classification_input_par_%(iter)d.par',
            'iter_par_block': 'Refine2D/Parameters/classification_input_par_%(iter)d_%(block)d.par',
            'iter_cls_block': 'Refine2D/ClassAverages/class_dump_file_%(iter)d_%(block)d.dump',
            'iter_cls_block_seed': 'Refine2D/ClassAverages/class_dump_file_%(iter)d_.dump'
        }
        self._updateFilenamesDict(myDict)

    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        parFn = self._getExtraPath(self._getFileName('iter_par',
                                                      iter=0))
        self._iterTemplate = parFn.replace('0', '*')
        self._iterRegex = re.compile('input_par_(\d{1,2})')

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
                      pointerCondition='hasCTF',
                      important=True, pointerClass='SetOfParticles',
                      help='Select the input particles.')
        form.addParam('inputClassAvg', PointerParam,
                      condition='not doContinue',
                      allowsNull=True,
                      label="Input class averages",
                      pointerClass='SetOfAverages',
                      help='Select starting class averages. If not provided, '
                           'they will be generated automatically.')
        form.addParam('areParticlesBlack', BooleanParam,
                       default=False,
                       label='Are the particles black?',
                       help='cisTEM requires particles to be black on white.')
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
        line.addParam('highResLimit1', FloatParam, default=40.0,
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
        form.addParam('angStep', FloatParam, default=15.0,
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
        line.addParam('rangeX', FloatParam, default=60.0,
                      label='X')
        line.addParam('rangeY', FloatParam, default=60.0,
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
        form.addParam('exclEdges', BooleanParam, default=False,
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
        form.addParam('percUsed', FloatParam, default=100.0,
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
        self._createFilenameTemplates()
        self._createIterTemplates()
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
        self._insertFunctionStep('convertInputStep')
        self.currPtcl = 1
        for iterN in self._allItersN():
            paramsDic = self._getParamsIteration(iterN)
            depsRefine = self._insertRefineIterStep(iterN, paramsDic)
            if iterN > 1:
                self._insertFunctionStep("mergeStep", iterN,
                                         prerequisites=depsRefine)

    def _insertRefineIterStep(self, iterN, paramsDic):
        """ Execute the refinement for the current iteration """
        depsRefine = []
        if iterN == 1:
            initParStepId = self._insertFunctionStep("writeInitParStep")
            refineId = self._insertFunctionStep("makeInitClassesStep",
                                                paramsDic,
                                                prerequisites=[initParStepId])
            depsRefine.append(refineId)
        else:
            jobs, ptcls_per_job = self._getJobsParams()
            for job in range(1, jobs + 1):
                refineId = self._insertFunctionStep("refineStep",
                                                    iterN, job,
                                                    ptcls_per_job,
                                                    paramsDic)
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

    def convertInputStep(self):
        self._createWorkingDirs()
        inputStack = self._getFileName('run_stack', run=0)
        inputClasses = self._getFileName('initial_cls')
        imgFn = self._getExtraPath(inputStack)
        self.writeParticlesByMic(imgFn)
        inputCls = self.inputClassAvg.get() if self.inputClassAvg else None

        if inputCls is not None:
            writeReferences(self.inputClassAvg.get(),
                            self._getExtraPath(inputClasses))

    def writeInitParStep(self):
        """ Construct a parameter file (.par). """
        #  This function will be called only for iterations 1 and 2.
        parFn = self._getExtraPath(self._getFileName('iter_par', iter=1))
        f = open(parFn, 'w')
        f.write("C           PSI   THETA     PHI       SHX       SHY     MAG  "
                "FILM      DF1      DF2  ANGAST  PSHIFT     OCC      LogP"
                "      SIGMA   SCORE  CHANGE\n")
        hasAlignment = self.hasAlignment()

        for i, part in self.iterParticlesByMic():
            ctf = part.getCTF()
            defocusU, defocusV, astig = ctf.getDefocusU(), ctf.getDefocusV(),\
                                        ctf.getDefocusAngle()
            phaseShift = ctf.getPhaseShift() or 0.00

            if hasAlignment:
                _, angles = geometryFromMatrix(part.getTransform().getMatrix())
                psi = angles[2]
            else:
                psi = 0.0

            string = '%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f'\
                     '%8.2f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n' % (
                i + 1, psi, 0., 0., 0., 0., 0, 0, defocusU, defocusV,
                astig, phaseShift, 100., 0, 10., 0., 0.)

            f.write(string)
        f.close()

    def makeInitClassesStep(self, paramsDic):
        argsStr = self._getRefineArgs()

        percUsed = self.numberOfClassAvg.get() * 300.0
        percUsed = percUsed / self._getPtclsNumber()  * 100.0
        if percUsed > 100.0:
            percUsed = 100.0

        paramsDic.update({
            'input_params': self._getFileName('iter_par', iter=1),
            'input_cls': '/dev/null',
            'output_cls': self._getFileName('iter_cls', iter=1),
            'output_params': '/dev/null',
            'percUsed': percUsed / 100.0,
            'dumpFn': '/dev/null'
        })

        cmdArgs = argsStr % paramsDic
        self.runJob(self._getProgram('refine2d'), cmdArgs,
                    cwd=self._getExtraPath(),
                    env=Plugin.getEnviron())

    def refineStep(self, iterN, job, ptcls_per_job, paramsDic):
        numPtcls = self._getPtclsNumber()

        if job == 1:
            firstPart = 1
            lastPart = 1 + int(ptcls_per_job)
            self.currPtcl = lastPart + 1
        else:
            firstPart = self.currPtcl
            lastPart = firstPart + int(ptcls_per_job)
            if lastPart > numPtcls:
                lastPart = numPtcls
            self.currPtcl = lastPart + 1

        argsStr = self._getRefineArgs()
        highRes = self._calcHighResLimit(self.finalIter,
                                         self.highResLimit1.get(),
                                         self.highResLimit2.get())

        percUsed = self._calcPercUsed(self.finalIter,
                                      iterN - 1,
                                      self.numberOfClassAvg.get(),
                                      numPtcls,
                                      self.percUsed.get(),
                                      self.autoPerc)

        paramsDic.update({
            'output_params': self._getFileName('iter_par_block', iter=iterN,
                                               block=job),
            'numberOfClassAvg': 0,  # determined from cls stack
            'firstPart': firstPart,
            'lastPart': lastPart,
            'percUsed': percUsed / 100.0,
            'highRes': highRes,
            'dump': 'YES',
            'dumpFn': self._getFileName('iter_cls_block', iter=iterN,
                                        block=job)
        })

        cmdArgs = argsStr % paramsDic
        self.runJob(self._getProgram('refine2d'), cmdArgs,
                    cwd=self._getExtraPath(),
                    env=Plugin.getEnviron())

    def mergeStep(self, iterN):
        jobs, _ = self._getJobsParams()
        self._mergeAllParFiles(iterN, jobs)
        argsStr = self._getMergeArgs()

        paramsDic = {
            'output_cls': self._getFileName('iter_cls', iter=iterN),
            'dumpSeed': self._getFileName('iter_cls_block_seed', iter=iterN),
            'numberOfJobs': jobs
        }

        cmdArgs = argsStr % paramsDic
        self.runJob(self._getProgram('merge2d'), cmdArgs,
                    cwd=self._getExtraPath(),
                    env=Plugin.getEnviron())

        dumpFns = self._getExtraPath('Refine2D/ClassAverages/class_dump_file_%d_*' % iterN)
        cleanPattern(dumpFns)

    def createOutputStep(self):
        partSet = self._getInputParticlesPointer()
        classes2D = self._createSetOfClasses2D(partSet.get())
        self._fillClassesFromIter(classes2D, self._lastIter())

        self._defineOutputs(outputClasses=classes2D)
        self._defineSourceRelation(partSet, classes2D)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if self.doContinue:
            errors.append('Continue mode not implemented yet!')

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
            self._getPtclsNumber())
        methods += "into %d classes using refine2d" % self.numberOfClassAvg
        return [methods]

    def _citations(self):
        return ['Sigworth1998', 'Scheres2005', 'Liang2015']

    # --------------------------- UTILS functions -----------------------------
    def _getProgram(self, program='refine2d'):
        return Plugin.getProgram(program)

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

    def _getPtclsNumber(self):
        return self._getInputParticles().getSize()

    def _getJobsParams(self):
        jobs = max(self.numberOfThreads.get(),
                   self.numberOfMpi.get())
        parts = self._getPtclsNumber()
        if parts - jobs < jobs:
            ptcls_per_job = 1.0
        else:
            ptcls_per_job = round(float(parts / jobs))

        return jobs, ptcls_per_job

    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1))  # group 1 is 1 digit iteration number
        return result

    def _lastIter(self):
        return self._getIterNumber(-1)

    def _fillClassesFromIter(self, clsSet, iterN):
        params = {'orderBy': ['_micId', 'id'],
                  'direction': 'ASC'
                  }
        self._classesInfo = {}  # store classes info, indexed by class id
        clsFn = self._getFileName('iter_cls', iter=iterN).replace('.mrc',
                                                                  '.mrc:mrcs')
        for classId in range(1, self.numberOfClassAvg.get() + 1):
            self._classesInfo[classId] = (classId,
                                          self._getExtraPath(clsFn))

        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=self._iterRows(iterN),
                             iterParams=params)

    def _updateParticle(self, item, row):
        vals = OrderedDict(zip(HEADER_COLUMNS, row))
        item.setClassId(vals.get('FILM'))
        item.setTransform(rowToAlignment(vals, item.getSamplingRate()))
        item._cistemLogP = em.Float(vals.get('LogP'))
        item._cistemSigma = em.Float(vals.get('SIGMA'))
        item._cistemOCC = em.Float(vals.get('OCC'))
        item._cistemScore = em.Float(vals.get('SCORE'))

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn = self._classesInfo[classId]
            item.getRepresentative().setLocation(index, fn)

    def _iterRows(self, iterN):
        filePar = self._getFileName('iter_par', iter=iterN)
        f1 = open(self._getExtraPath(filePar))
        for line in f1:
            if not line.startswith("C"):
                values = map(float, line.strip().split())

                yield values

        f1.close()

    def iterParticlesByMic(self):
        """ Iterate the particles ordered by micrograph """
        for i, part in enumerate(self._getInputParticles().iterItems(orderBy=['_micId', 'id'],
                                                                     direction='ASC')):
            yield i, part

    def writeParticlesByMic(self, stackFn):
        self._getInputParticles().writeStack(stackFn,
                                             orderBy=['_micId', 'id'],
                                             direction='ASC')

    def hasAlignment(self):
        inputParts = self._getInputParticles()
        return inputParts.hasAlignment()

    def _getParamsIteration(self, iterN):
        """ Defining the current iteration """
        imgSet = self._getInputParticles()
        acq = imgSet.getAcquisition()

        # Prepare arguments to call refine2d
        paramsDic = {'input_stack': self._getFileName('run_stack', run=0),
                     'input_params': self._getFileName('iter_par', iter=iterN-1),
                     'input_cls': self._getFileName('iter_cls', iter=iterN-1),
                     'output_params': self._getFileName('iter_par', iter=iterN),
                     'output_cls': self._getFileName('iter_cls', iter=iterN),
                     'numberOfClassAvg': self.numberOfClassAvg.get(),
                     'firstPart': 1,
                     'lastPart': 0,
                     'percUsed': self.percUsed.get() / 100.0,
                     'pixSize': imgSet.getSamplingRate(),
                     'voltage': acq.getVoltage(),
                     'sphAber': acq.getSphericalAberration(),
                     'ampCont': acq.getAmplitudeContrast(),
                     'maskRad': self.maskRad.get(),
                     'lowRes': self.lowResLimit.get(),
                     'highRes': self.highResLimit1.get(),
                     'angStep': self.angStep.get(),
                     'rangeX': self.rangeX.get(),
                     'rangeY': self.rangeY.get(),
                     'smooth': self.smooth.get(),
                     'pad': 2,
                     'normalize': 'YES',
                     'invertContrast': 'NO' if self.areParticlesBlack else 'YES',
                     'exclEdges': 'YES' if self.exclEdges else 'NO',
                     'dump': 'NO',
                     'dumpFn': self._getFileName('iter_cls_block', iter=iterN,
                                                 block=1)
                     }

        return paramsDic

    def _getRefineArgs(self):
        argsStr = """ << eof
%(input_stack)s
%(input_params)s
%(input_cls)s
%(output_params)s
%(output_cls)s
%(numberOfClassAvg)d
%(firstPart)d
%(lastPart)d
%(percUsed)f
%(pixSize)f
%(voltage)f
%(sphAber)f
%(ampCont)f
%(maskRad)f
%(lowRes)f
%(highRes)f
%(angStep)f
%(rangeX)f
%(rangeY)f
%(smooth)f
%(pad)d
%(normalize)s
%(invertContrast)s
%(exclEdges)s
%(dump)s
%(dumpFn)s
eof
"""
        return argsStr

    def _mergeAllParFiles(self, iterN, numberOfBlocks):
        """ This method merge all parameters files that has been
        created in a refineStep. """
        self._enterDir(self._getExtraPath())
        outFn = self._getFileName('iter_par', iter=iterN)

        if numberOfBlocks != 1:
            f1 = open(outFn, 'w+')
            f1.write("C           PSI   THETA     PHI       SHX       SHY     MAG  "
                "FILM      DF1      DF2  ANGAST  PSHIFT     OCC      LogP"
                "      SIGMA   SCORE  CHANGE\n")
            for block in range(1, numberOfBlocks + 1):
                parFn = self._getFileName('iter_par_block', iter=iterN,
                                          block=block)
                if not exists(parFn):
                    raise Exception("Error: file %s does not exist" % parFn)
                f2 = open(parFn)

                for l in f2:
                    if not l.startswith('C'):
                        f1.write(l)
                f2.close()
                cleanPattern(parFn)
            f1.close()
        else:
            parFn = self._getFileName('iter_par_block', iter=iterN, block=1)
            moveFile(parFn, outFn)

        self._leaveDir()

    def _getMergeArgs(self):
        argsStr = """ << eof
%(output_cls)s
%(dumpSeed)s
%(numberOfJobs)d
eof
"""
        return argsStr

    def _calcHighResLimit(self, iter_total, highRes1, highRes2):
        """ Ramp up the high resolution limit.
        Copied from MyRefine2DPanel.cpp of cisTEM. """
        if iter_total > 1:
            if iter_total >= 4:
                reachMaxHighResAtCycle = iter_total * 3 / 4
            else:
                reachMaxHighResAtCycle = iter_total
            if iter_total >= reachMaxHighResAtCycle:
                highRes = highRes2
            else:
                highRes = highRes1 + iter_total / (reachMaxHighResAtCycle - 1) * (highRes2 - highRes1)
        else:
            highRes = highRes2

        return highRes

    def _calcPercUsed(self, iter_total, iter_done, numCls,
                      numParts, percUsed, autoPerc=True):
        """ Copied from MyRefine2DPanel.cpp of cisTEM. """
        minPercUsed = percUsed
        if autoPerc:
            if iter_total < 10:
                percUsed = 100.0
            else:
                if iter_total < 20:
                    if iter_done < 5:
                        percUsed = numCls * 300 / numParts * 100.0
                        if percUsed > 100.0:
                            percUsed = 100.0
                    else:
                        if iter_done < iter_total - 5:
                            percUsed = numCls * 300 / numParts * 100.0
                            if percUsed > 100.0:
                                percUsed = 100.0
                            elif percUsed < 30:
                                percUsed = 30.0
                        else:
                            percUsed = 100.0
                elif iter_total < 30:
                    if iter_done < 10:
                        percUsed = numCls * 300 / numParts * 100.0
                        if percUsed > 100.0:
                            percUsed = 100.0
                    elif iter_done < iter_total - 5:
                        percUsed = numCls * 300 / numParts * 100.0
                        if percUsed > 100.0:
                            percUsed = 100.0
                        elif percUsed < 30:
                            percUsed = 30.0
                    else:
                        percUsed = 100.0
                else:
                    if iter_done < 15:
                        percUsed = numCls * 300 / numParts * 100.0
                        if percUsed > 100.0:
                            percUsed = 100.0
                    elif iter_done < iter_total - 5:
                        percUsed = numCls * 300 / numParts * 100.0
                        if percUsed > 100.0:
                            percUsed = 100.0
                        elif percUsed < 30:
                            percUsed = 30.0
                    else:
                        percUsed = 100.0
        else:
            percUsed = percUsed
        if percUsed < minPercUsed:
            percUsed = minPercUsed

        return percUsed
