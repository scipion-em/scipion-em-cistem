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

from pyworkflow.protocol import STEPS_PARALLEL
import pyworkflow.em as em
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        IntParam, StringParam,
                                        BooleanParam)
import cistem
from cistem.convert import *



class CistemProtAbInitio(em.ProtInitialVolume):
    """ Protocol to run ab-initio reconstruction in cisTEM. """
    _label = 'ab-initio'

    def __init__(self, **args):
        em.ProtInitialVolume.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasCTF',
                      important=True,
                      label="Input particles",
                      help='Select the input particles.')
        form.addParam('numOfStarts', IntParam,
                      default=2,
                      label='Number of starts',
                      help='The number of times the ab-initio reconstruction '
                           'is restarted, using the result from the previous '
                           'run in each restart.')
        form.addParam('numOfCycles', IntParam,
                      default=40,
                      label='Number of cycles per start',
                      help='The number of refinement cycles to run for '
                           'each start. The percentage of particles and '
                           'the refinement resolution limit will be '
                           'adjusted automatically from cycle to cycle '
                           'using initial and final values specified '
                           'under Expert Options.')

        form.addSection(label='Expert options')
        group = form.addGroup('Refinement')
        group.addParam('StartResLimit', FloatParam, default=20.0,
                      label='Initial resolution limit (A)',
                      help='The starting resolution limit used to align '
                           'particles against the current 3D '
                           'reconstruction. In most cases, this should '
                           'specify a relatively low resolution to make '
                           'sure the reconstructions generated in the '
                           'initial refinement cycles do not develop '
                           'spurious high-resolution features.')
        group.addParam('EndResLimit', FloatParam, default=8.0,
                       label='Final resolution limit (A)',
                       help='The resolution limit used in the final '
                            'refinement cycle. In most cases, this '
                            'should specify a resolution at which '
                            'expected secondary structure becomes '
                            'apparent, i.e. around 9 A.')

        line = group.addLine('Mask radius (A):',
                             help='The radius of the circular mask applied. '
                           'This mask should be sufficiently large to include '
                           'the largest dimension of the particle. The mask '
                           'helps remove noise outside the area of the '
                           'particle.')
        line.addParam('maskRadInner', FloatParam, default=0.0,
                      label='inner')
        line.addParam('maskRad', FloatParam, default=90.0,
                      label='outer')

        line = group.addLine('Search range (A): ',
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

        group.addParam('autoMask', BooleanParam, default=True,
                      label='Use auto-masking?',
                      help='Should the 3D reconstructions be masked? '
                           'Masking is important to suppress weak '
                           'density features that usually appear in '
                           'the early stages of ab-initio reconstruction, '
                           'thus preventing them to get amplified during '
                           'the iterative refinement. Masking should only '
                           'be disabled if it appears to interfere with '
                           'the reconstruction process.')
        group.addParam('autoPerc', BooleanParam, default=True,
                      label='Auto percent used?',
                      help='Should the percentage of particles used in '
                           'each refinement cycle be set automatically? '
                           'If reconstructions appear very noisy or '
                           'reconstructions settle into a wrong structure '
                           'that does not change anymore during iterations, '
                           'disable this option and specify initial and final '
                           'percentages manually. To reduce noise, increase '
                           'the percentage; to make reconstructions more '
                           'variable, decrease the percentage. By default, '
                           'the initial percentage is set to include an '
                           'equivalent of 2500 asymmetric units and the '
                           'final percentage corresponds to 10,000 asymmetric '
                           'units used.')

        line = group.addLine('Percent used (%):',
                             condition='not autoPerc',
                             help='User-specified percentages of particles '
                                  'used when Auto Percent Used is disabled.')
        line.addParam('percUsed1', FloatParam, default=10.0,
                      condition='not autoPerc',
                      label='initial')
        line.addParam('percUsed2', FloatParam, default=10.0,
                      condition='not autoPerc',
                      label='final')

        group.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry",
                      help='')
        group.addParam('applySym', BooleanParam, default=False,
                       label='Always apply symmetry?',
                       help='')

        group = form.addGroup('Reconstruction')
        group.addParam('applyBlur', BooleanParam, default=False,
                       label='Apply likelihood blurring?',
                       help='Should the reconstructions be blurred by '
                            'inserting each particle image at multiple '
                            'orientations, weighted by a likelihood '
                            'function? Enable this option if the ab-initio '
                            'procedure appears to suffer from over-fitting '
                            'and the appearance of spurious high-resolution '
                            'features.')
        group.addParam('smooth', FloatParam, default=1.0,
                      condition='applyBlur',
                      label='Smoothing factor [0-1]',
                      help='A factor that reduces the range of likelihoods '
                           'used for blurring. A smaller number leads to more '
                           'blurring. The user should try values between '
                           '0.1 and 1.')

        form.addParallelSection(threads=4, mpi=1)

    # -------------------------- INSERT steps functions -----------------------

    # -------------------------- STEPS functions ------------------------------

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        summary = []

        return summary

    def _citations(self):
        return ['Grigorieff2016']

    # -------------------------- UTILS functions ------------------------------
