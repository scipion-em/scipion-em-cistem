# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
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

import pyworkflow.protocol.params as params
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.em.constants import RELATION_CTF
from pyworkflow.em import ProtParticlePickingAuto
from pyworkflow.em.convert import ImageHandler
import pyworkflow.utils as pwutils
from pyworkflow.utils.properties import Message

from cistem import Plugin
from cistem.convert import readSetOfCoordinates, writeReferences
from cistem.constants import *


class CistemProtFindParticles(ProtParticlePickingAuto):
    """ Protocol to pick particles in a set of micrographs using cisTEM. """
    _label = 'find particles'

    def __init__(self, **kwargs):
        ProtParticlePickingAuto.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        ProtParticlePickingAuto._defineParams(self, form)
        form.addParam('ctfRelations', params.RelationParam,
                      relationName=RELATION_CTF,
                      attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to the '
                           'input micrographs.')
        form.addParam('pickType', params.EnumParam, default=0,
                      important=True,
                      label='Picking algorithm',
                      choices=['Ab-initio', 'Reference-based'],
                      display=params.EnumParam.DISPLAY_HLIST)
        form.addParam('inputRefs', params.PointerParam,
                      condition='pickType==1', important=True,
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      label='Input references',
                      help='Provide a set of 2D templates to use in '
                           'the search.')
        form.addParam('maxradius', params.FloatParam, default=120.0,
                      label='Max particle radius (A)',
                      help='In Angstroms, the maximum radius of the '
                            'particles to be found. This also determines '
                            'the minimum distance between picks.')
        form.addParam('radius', params.FloatParam, default=80.0,
                      label='Characteristic particle radius (A)',
                      help='In Angstroms, the radius within which most '
                           'of the density is enclosed. The template '
                           'for picking is a soft-edge disc, where '
                           'the edge is 5 pixels wide and this '
                           'parameter defines the radius at which '
                           'the cosine-edge template reaches 0.5.')
        form.addParam('threshold', params.FloatParam, default=6.0,
                      label='Threshold peak height',
                      help='Particle coordinates will be defined as '
                           'the coordinates of any peak in the search '
                           'function which exceeds this threshold. '
                           'In numbers of standard deviations '
                           'above expected noise variations in '
                           'the scoring function. See Sigworth '
                           '(2004) for definition.')
        form.addParam('avoidHighVar', params.BooleanParam, default=False,
                      label='Avoid high variance areas',
                      help='Avoid areas with abnormally high local variance. '
                           'This can be effective in avoiding edges '
                           'of support films or contamination.')
        form.addParam('ptclWhite', params.BooleanParam, default=False,
                      label='Particles are white on a dark background?')

        form.addSection(label='Expert options')
        form.addParam('highRes', params.FloatParam, default=30.0,
                      label='Highest resolution used in picking (A)',
                      help='The template and micrograph will be resampled '
                           '(by Fourier cropping) to a pixel size of half '
                           'the resolution given here. Note that the '
                           'information in the corners of the Fourier '
                           'transforms (beyond the Nyquist frequency) '
                           'remains intact, so that there is some small '
                           'risk of bias beyond this resolution.')
        form.addParam('minDist', params.IntParam, default=0,
                      label='Minimum distance from edges (px)',
                      help='No particle shall be picked closer than '
                           'this distance from the edges of the micrograph. '
                           'In pixels.')
        form.addParam('useRadAvg', params.BooleanParam, default=True,
                      condition='pickType==1',
                      label='Use radial averages of templates',
                      help='Say yes if the templates should be '
                           'rotationally averaged')
        form.addParam('rotateRef', params.IntParam, default=0,
                      condition='pickType==1 and not useRadAvg',
                      label='Rotate each template this many times',
                      help='If > 0, each template image will be '
                           'rotated this number of times and the '
                           'micrograph will be searched for the rotated '
                           'template.')
        form.addParam('avoidLocMean', params.BooleanParam, default=True,
                      label='Avoid areas with abnormal local mean',
                      help='Avoid areas with abnormally low or high '
                           'local mean. This can be effective to avoid '
                           'picking from, e.g., contaminating ice crystals, '
                           'support film.')
        form.addParam('bgBoxes', params.IntParam, default=30,
                      label='Number of background boxes',
                      help='Number of background areas to use in estimating '
                           'the background spectrum. The larger the number '
                           'of boxes, the more accurate the estimate should '
                           'be, provided that none of the background boxes '
                           'contain any particles to be picked.')
        form.addParam('bgAlgo', params.EnumParam, default=LOW_VARIANCE,
                      choices=['Lowest variance', 'Variance near mode'],
                      display=params.EnumParam.DISPLAY_COMBO,
                      label='Algorithm to find background areas',
                      help='Testing so far suggests that areas of lowest '
                           'variance in experimental micrographs should be '
                           'used to estimate the background spectrum. '
                           'However, when using synthetic micrographs this '
                           'can lead to bias in the spectrum estimation '
                           'and the alternative (areas with local variances '
                           'near the mean of the distribution of local '
                           'variances) seems to perform better')

        self._defineStreamingParams(form)

        # Allow many threads if we can put more than one in a CPU
        form.addParallelSection(threads=1, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self.inputStreaming = self.getInputMicrographs().isStreamOpen()

        if self.streamingBatchSize > 0 or self.inputStreaming:
            # If the input is in streaming, follow the base class policy
            # about inserting new steps and discovery new input/output
            ProtParticlePickingAuto._insertAllSteps(self)
            self.createOutputStep = self._doNothing
        else:
            # If not in streaming, then we will just insert a single step to
            # pick all micrographs at once since it is much faster
            self._insertInitialSteps()
            self._insertFunctionStep('_pickMicrographStep',
                                     self.getInputMicrographs(),
                                     *self._getPickArgs())
            self._insertFunctionStep('createOutputStep')

            # Disable streaming functions:
            self._insertFinalSteps = self._doNothing
            self._stepsCheck = self._doNothing

    def _insertInitialSteps(self):
        # Convert the input micrographs to mrc
        inputRefs = self.getInputReferences()
        refsId = inputRefs.strId() if inputRefs is not None else None
        convertId = self._insertFunctionStep('convertInputStep',
                                             self.getInputMicrographs().strId(),
                                             refsId)
        return [convertId]

    def _doNothing(self, *args):
        pass  # used to avoid some streaming functions

    def _loadInputList(self):
        """ This function is re-implemented in this protocol, because it has
         a SetOfCTF as input, so for streaming, we only want to report the
         micrographs for which the CTF is ready.
        """
        micDict, micClose = self._loadMics(self.getInputMicrographs())
        ctfDict, ctfClosed = self._loadCTFs(self.ctfRelations.get())

        # Remove the micrographs that have no CTF
        # and set the CTF property for those who have it
        for micKey, mic in micDict.iteritems():
            if micKey in ctfDict:
                mic.setCTF(ctfDict[micKey])
            else:
                del micDict[micKey]

        # Return the updated micDict and the closed status
        return micDict, micClose and ctfClosed

    # --------------------------- STEPS functions -------------------------------
    def convertInputStep(self, micsId, refsId):
        # Match ctf information against the micrographs
        self.ctfDict = {}
        if self.ctfRelations.get() is not None:
            for ctf in self.ctfRelations.get():
                self.ctfDict[ctf.getMicrograph().getMicName()] = ctf.clone()

        ih = ImageHandler()
        for mic in self.getInputMicrographs():
            micName = mic.getFileName()
            # We convert the input micrographs if they are not .mrc
            outMic = os.path.join(self._getTmpPath(),
                                  pwutils.replaceBaseExt(micName, 'mrc'))
            if micName.endswith('.mrc'):
                pwutils.createLink(micName, outMic)
            else:
                ih.convert(micName, outMic)

        if refsId is not None:
            writeReferences(self.getInputReferences(),
                            self._getExtraPath('references.mrc'))

    def _pickMicrograph(self, mic, args):
        self._pickMicrographStep([mic], args)

    def _pickMicrographList(self, micList, args):
        self._pickMicrographStep(micList, args)

    def _pickMicrographStep(self, mics, args):
        for mic in mics:
            micName = mic.getFileName()
            outMic = os.path.join(self._getTmpPath(),
                                  pwutils.replaceBaseExt(micName, 'mrc'))
            ctf = self.ctfDict[mic.getMicName()]

            args.update({'micName': outMic,
                         'logFn': self._getLogFn(mic),
                         'outStack': self._getStackFn(mic),
                         'phaseShift': ctf.getPhaseShift() or 0.0,
                         'defocusU': ctf.getDefocusU(),
                         'defocusV': ctf.getDefocusV(),
                         'defocusAngle': ctf.getDefocusAngle()
                         })

            if self.pickType == 1:
                args.update({
                    'refsFn': self._getExtraPath('references.mrc'),
                    'useRadAvg': 'YES' if self.useRadAvg else 'NO',
                    'rotateRef': self.rotateRef.get(),
                })

            argsStr = self._getArgsStr()
            cmdArgs = argsStr % args

            try:
                self.runJob(self._getProgram(), cmdArgs,
                            env=Plugin.getEnviron())

                # Move output from micPath (tmp) to extra
                pltFn = pwutils.replaceExt(self._getStackFn(mic), 'plt')
                pwutils.moveFile(pltFn, self._getPltFn(mic))

                # Clean tmp folder
                pwutils.cleanPath(outMic)
                pwutils.cleanPath(self._getLogFn(mic))
            except:
                print("Picking for mic %s failed! Skipping.." % micName)

    def createOutputStep(self):
        micSet = self.getInputMicrographs()
        coordSet = self._createSetOfCoordinates(micSet)
        self.readCoordsFromMics(self._getExtraPath(), micSet,
                                coordSet)

        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(self.inputMicrographs, coordSet)

    # --------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        summary.append("Number of input micrographs: %d"
                       % self.getInputMicrographs().getSize())
        if self.getOutputsSize() > 0:
            summary.append("Number of particles picked: %d"
                           % self.getCoords().getSize())
            summary.append("Particle size: %d" % self.getCoords().getBoxSize())
            summary.append("Threshold: %0.2f" % self.threshold)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def _citations(self):
        return ['Sigworth2004']

    def _methods(self):
        methodsMsgs = []
        if self.getInputMicrographs() is None:
            return ['Input micrographs not available yet.']
        methodsMsgs.append("Input micrographs %s of size %d."
                           % (self.getObjectTag(self.getInputMicrographs()),
                              self.getInputMicrographs().getSize()))

        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append('%s: User picked %d particles with a particle '
                               'size of %d and threshold %0.2f.'
                               % (self.getObjectTag(output), output.getSize(),
                                  output.getBoxSize(), self.threshold))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    # --------------------------- UTILS functions -------------------------------
    def _getProgram(self):
        return Plugin.getProgram(FIND_PARTICLES_BIN)

    def _getPickArgs(self):
        """ Format arguments to call find_particles program. """
        inputMics = self.getInputMicrographs()
        acq = inputMics.getAcquisition()
        sampling = inputMics.getSamplingRate()

        self.argsDict = {'samplingRate': sampling,
                         'voltage': acq.getVoltage(),
                         'cs': acq.getSphericalAberration(),
                         'ampContrast': acq.getAmplitudeContrast(),
                         'templates': 'NO' if self.pickType == 0 else 'YES',
                         'radius': self.radius.get(),
                         'maxradius': self.maxradius.get(),
                         'highRes': self.highRes.get(),
                         'boxSize': 120,
                         'minDist': self.minDist.get(),
                         'threshold': self.threshold.get(),
                         'avoidHighVar': 'YES' if self.avoidHighVar else 'NO',
                         'avoidLocMean': 'YES' if self.avoidLocMean else 'NO',
                         'algorithm': self.bgAlgo.get(),
                         'bgBoxes': self.bgBoxes.get(),
                         'ptclWhite': 'YES' if self.ptclWhite else 'NO'
                         }

        return [self.argsDict]

    def _getArgsStr(self):
        argsStr = """ << eof > %(logFn)s
%(micName)s
%(samplingRate)f
%(voltage)f
%(cs)f
%(ampContrast)f
%(phaseShift)f
%(defocusU)f
%(defocusV)f
%(defocusAngle)f"""
        if self.pickType == 0:
            argsStr += """
NO
%(radius)f
%(maxradius)f
%(highRes)f
%(outStack)s
%(boxSize)d
%(minDist)d
%(threshold)f
%(avoidHighVar)s
%(avoidLocMean)s
%(algorithm)d
%(bgBoxes)d
%(ptclWhite)s
eof"""
        else: # ref-based picking
            argsStr += """
YES
%(refsFn)s
%(useRadAvg)s"""

            if self.useRadAvg:
                argsStr += """
%(maxradius)f
%(highRes)f
%(outStack)s
%(boxSize)d
%(minDist)d
%(threshold)f
%(avoidHighVar)s
%(avoidLocMean)s
%(algorithm)d
%(bgBoxes)d
%(ptclWhite)s
eof"""
            else:
                if self.rotateRef > 0:
                    argsStr += """
YES
%(rotateRef)d
%(maxradius)f
%(highRes)f
%(outStack)s
%(boxSize)d
%(minDist)d
%(threshold)f
%(avoidHighVar)s
%(avoidLocMean)s
%(algorithm)d
%(bgBoxes)d
%(ptclWhite)s
eof"""
                else:
                    argsStr += """
NO
%(maxradius)f
%(highRes)f
%(outStack)s
%(boxSize)d
%(minDist)d
%(threshold)f
%(avoidHighVar)s
%(avoidLocMean)s
%(algorithm)d
%(bgBoxes)d
%(ptclWhite)s
eof"""

        return argsStr

    def readCoordsFromMics(self, extraDir, micList, coordSet):
        if coordSet.getBoxSize() is None:
            coordSet.setBoxSize(self._getBoxSize())

        readSetOfCoordinates(self._getExtraPath(), micList, coordSet)

    def _getBoxSize(self):
        return 128

    def _getMicrographDir(self, mic):
        """ Return an unique dir name for results of the micrograph. """
        return self._getTmpPath('mic_%06d' % mic.getObjId())

    def getInputMicrographs(self):
        return self.inputMicrographs.get()

    def _getLogFn(self, mic):
        micName = mic.getFileName()
        return pwutils.join(self._getTmpPath(),
                            pwutils.replaceBaseExt(micName, 'log'))

    def _getStackFn(self, mic):
        return self._getTmpPath('mic_%06d.mrc' % mic.getObjId())

    def _getPltFn(self, mic):
        micName = mic.getFileName()
        return pwutils.join(self._getExtraPath(),
                            pwutils.replaceBaseExt(micName, 'plt'))

    def getInputReferences(self):
        return self.inputRefs.get() if self.inputRefs else None
