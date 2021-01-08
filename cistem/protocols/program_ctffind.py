# **************************************************************************
# *
# * Authors:     Josue Gomez BLanco (josue.gomez-blanco@mcgill.ca) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] Department of Anatomy and Cell Biology, McGill University
# * [2] SciLifeLab, Stockholm University
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

from numpy import deg2rad

from pwem.objects import CTFModel
import pyworkflow.protocol.params as params

from cistem import Plugin
from ..constants import CTFFIND4_BIN
from ..convert import parseCtffind4Output, readCtfModel


class ProgramCtffind:
    """
    Wrapper of Ctffind4 program that will handle parameters definition
    and also execution of the program with the proper arguments.
    This class is not a Protocol, but somewhat related, since it can be used from
    protocols that perform CTF estimation.
    """
    def __init__(self, protocol):
        self._program = Plugin.getProgram(CTFFIND4_BIN)  # Load program to use
        self._findPhaseShift = protocol.findPhaseShift.get()
        self._args, self._params = self._getArgs(protocol)  # Load general arguments

    @classmethod
    def defineInputParams(cls, form):
        """ Define input/common parameters. """
        form.addSection(label='Input')
        form.addParam('recalculate', params.BooleanParam, default=False,
                      condition='recalculate',
                      label="Do recalculate ctf?")
        form.addParam('continueRun', params.PointerParam, allowsNull=True,
                      condition='recalculate', label="Input previous run",
                      pointerClass='CistemProtCTFFind')
        form.addHidden('sqliteFile', params.FileParam, condition='recalculate',
                       allowsNull=True)

        form.addParam('inputType', params.EnumParam, default=1,
                      label='Estimate using:',
                      choices=['Movies', 'Micrographs'],
                      display=params.EnumParam.DISPLAY_HLIST)
        form.addParam('inputMicrographs', params.PointerParam, important=True,
                      condition='not recalculate and inputType==1',
                      label='Input micrographs',
                      pointerClass='SetOfMicrographs')
        form.addParam('inputMovies', params.PointerParam, important=True,
                      condition='not recalculate and inputType==0',
                      label='Input movies',
                      pointerClass='SetOfMovies')
        form.addParam('avgFrames', params.IntParam, default=3,
                      condition='inputType==0',
                      label='No. movie frames to average',
                      help='When estimating parameters from movie frames, '
                           'enter how many frames should be included '
                           'in the sub-averages used to calculate '
                           'the amplitude spectra.')
        form.addParam('usePowerSpectra', params.BooleanParam, default=False,
                      condition='not recalculate and inputType==1',
                      label="Use power spectra?",
                      help="If set to Yes, the CTF estimation will be done "
                           "using power spectra calculated during "
                           "Relion motion correction.")

    @classmethod
    def defineProcessParams(cls, form):
        """ Define specific parameters. """
        form.addParam('windowSize', params.IntParam, default=512,
                      label='FFT box size (px)', condition='not recalculate',
                      help='The dimensions (in pixels) of the amplitude '
                           'spectrum CTFfind will compute. Smaller box '
                           'sizes make the fitting process significantly '
                           'faster, but sometimes at the expense of '
                           'fitting accuracy. If you see warnings '
                           'regarding CTF aliasing, consider '
                           'increasing this parameter.')

        group = form.addGroup('Search limits')
        line = group.addLine('Resolution (A)', condition='not recalculate',
                             help='The CTF model will be fit to regions '
                                  'of the amplitude spectrum corresponding '
                                  'to this range of resolution.')
        line.addParam('lowRes', params.FloatParam, default=30., label='Min')
        line.addParam('highRes', params.FloatParam, default=5., label='Max')

        line = group.addLine('Defocus search range (A)',
                             condition='not recalculate',
                             help='Select _minimum_ and _maximum_ values for '
                                  'defocus search range (in A). Underfocus'
                                  ' is represented by a positive number.')
        line.addParam('minDefocus', params.FloatParam, default=5000.,
                      label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=50000.,
                      label='Max')
        group.addParam('stepDefocus', params.FloatParam, default=500.,
                       label='Defocus step (A)',
                       help='Step size for the defocus search.')

        group.addParam('slowSearch', params.BooleanParam, default=False,
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Slower, more exhaustive search?",
                       help="Select this option if CTF determination "
                            "fails on images that show clear Thon rings "
                            "and should therefore yield good CTF parameters, "
                            "or if you expect noticeably elliptical Thon "
                            "rings and high noise.")

        group.addParam('fixAstig', params.BooleanParam, default=True,
                       label='Restrain astigmatism?',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Should the amount of astigmatism be restrained '
                            'during the parameter search and refinement? '
                            'This option should be selected when astigmatism '
                            'is expected to be small to produce more reliable '
                            'fits. Disable this option if you expect '
                            'large astigmatism.')
        group.addParam('astigmatism', params.FloatParam,
                       default=100.0, condition='fixAstig',
                       label='Tolerated astigmatism (A)',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='When restraining astigmatism, astigmatism values '
                            'much larger than this will be penalized. '
                            'Set to negative to remove this restraint. '
                            'In cases where the amplitude spectrum is '
                            'very noisy, such a restraint can help '
                            'achieve more accurate results.')

        form.addSection(label='Phase shift')
        form.addParam('findPhaseShift', params.BooleanParam,
                      default=False,
                      label="Find additional phase shift?",
                      help='Was the data recorded using a phase plate '
                           'with variable phase shift that must be '
                           'determined together with the defocus '
                           'parameters?')
        form.addParam('minPhaseShift', params.FloatParam, default=0.,
                      label="Minimum phase shift (deg)",
                      condition='findPhaseShift',
                      help='If finding an additional phase shift, '
                           'this value sets the lower bound for the '
                           'search.')
        form.addParam('maxPhaseShift', params.FloatParam, default=180,
                      label="Maximum phase shift (deg)",
                      condition='findPhaseShift',
                      help='If finding an additional phase shift, '
                           'this value sets the upper bound for the '
                           'search.')
        form.addParam('stepPhaseShift', params.FloatParam, default=10,
                      label="Phase shift search step (deg)",
                      condition='findPhaseShift',
                      help='If finding an additional phase shift, '
                           'this value sets the step size for the '
                           'search.')

        form.addParallelSection(threads=2, mpi=1)

    def getCommand(self, **kwargs):
        """
        :param kwargs: The input keywords argument should contain key-values
        for one micrograph or group of micrographs.
        :return: the program and arguments to be run
        """
        paramDict = dict(self._params)
        paramDict.update(kwargs)
        return self._program, self._args % paramDict

    def parseOutput(self, filename):
        """ Retrieve defocus U, V and angle from the
        output file of the program execution.
        :param filename: input file to parse
        :return: a tuple of CTF values
        """
        return parseCtffind4Output(filename)

    def parseOutputAsCtf(self, filename, psdFile=None):
        """ Parse the output file and build the CTFModel object
        with the values.
        :param filename: input file to parse
        :param psdFile: if defined, set PSD for the CTF model
        """
        ctf = CTFModel()
        readCtfModel(ctf, filename)
        if psdFile:
            ctf.setPsdFile(psdFile)

        return ctf

    def _getArgs(self, protocol):
        """ Update first the params dict.
        :param protocol: input protocol instance
        :return: args string and matching params dict
        """
        paramDict = protocol.getCtfParamsDict()
        paramDict['step_focus'] = protocol.stepDefocus.get()
        paramDict['fixAstig'] = "yes" if protocol.fixAstig else "no"
        paramDict['astigmatism'] = protocol.astigmatism.get()
        paramDict['lowRes'] = protocol.lowRes.get()
        if paramDict['lowRes'] > 50:
            paramDict['lowRes'] = 50
        paramDict['highRes'] = protocol.highRes.get()
        # defocus is in Angstroms now
        paramDict['minDefocus'] = protocol.minDefocus.get()
        paramDict['maxDefocus'] = protocol.maxDefocus.get()

        if self._findPhaseShift:
            paramDict['phaseShift'] = "yes"
            paramDict['minPhaseShift'] = deg2rad(protocol.minPhaseShift.get())
            paramDict['maxPhaseShift'] = deg2rad(protocol.maxPhaseShift.get())
            paramDict['stepPhaseShift'] = deg2rad(protocol.stepPhaseShift.get())
        else:
            paramDict['phaseShift'] = "no"

        paramDict['slowSearch'] = "yes" if protocol.slowSearch else "no"

        args = """   << eof > %(ctffindOut)s
%(micFn)s
%(ctffindPSD)s
%(samplingRate)f
%(voltage)f
%(sphericalAberration)f
%(ampContrast)f
%(windowSize)d
%(lowRes)f
%(highRes)f
%(minDefocus)f
%(maxDefocus)f
%(step_focus)f
no
%(slowSearch)s
%(fixAstig)s
%(phaseShift)s
no
eof\n
"""

        if protocol.usePowerSpectra:
            args = args.replace('<< eof > %(ctffindOut)s',
                                '--amplitude-spectrum-input << eof > %(ctffindOut)s')
            args = args.replace('%(samplingRate)f',
                                '%(powerSpectraPix)f')

        if protocol.fixAstig:
            args = args.replace('%(fixAstig)s',
                                '%(fixAstig)s\n'
                                '%(astigmatism)f')

        if self._findPhaseShift:
            args = args.replace('%(phaseShift)s',
                                '%(phaseShift)s\n'
                                '%(minPhaseShift)f\n'
                                '%(maxPhaseShift)f\n'
                                '%(stepPhaseShift)f')
        return args, paramDict
