# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca) [1]
# *
# * [1] Department of Anatomy and Cell Biology, McGill University
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
import numpy as np
from collections import OrderedDict
import logging
logger = logging.getLogger(__name__)

from pwem.objects import (Coordinate, SetOfClasses2D, SetOfAverages,
                          Transform, CTFModel, Float)
from pwem.constants import ALIGN_PROJ
from pwem.emlib.image import ImageHandler
import pwem.convert.transformations as transformations
from pyworkflow.utils import replaceBaseExt, exists


HEADER_COLUMNS = ['INDEX', 'PSI', 'THETA', 'PHI', 'SHX', 'SHY', 'MAG',
                  'FILM', 'DF1', 'DF2', 'ANGAST', 'PSHIFT', 'OCC',
                  'LogP', 'SIGMA', 'SCORE', 'CHANGE']


class FrealignParFile(object):
    """ Handler class to read/write Frealign par file."""
    def __init__(self, filename, mode='r'):
        self._file = open(filename, mode)
        self._count = 0

    def __iter__(self):
        """ Convert a line into a dict with HEADER_COLUMNS as keys.
        :return: yield a dict - single row
        """
        for line in self._file:
            line = line.strip()
            if not line.startswith('C'):
                row = OrderedDict(zip(HEADER_COLUMNS, line.split()))
                yield row

    def close(self):
        self._file.close()


def readSetOfParticles(inputSet, outputSet, parFileName):
    """ Iterate through the inputSet and the parFile lines
     and populate the outputSet with the same particles
     of inputSet, but with the angles and shift (3d alignment)
     updated from the parFile info.
     It is assumed that the order of iteration of the particles
     and the lines match and have the same number.
     :param inputSet: input set of particles
     :param outputSet: output set of particles to be populated
     :param parFileName: Frealign par file to read alignments
     """
    # create dictionary that matches input particles with param file
    samplingRate = inputSet.getSamplingRate()
    parFile = FrealignParFile(parFileName)
    partIter = iter(inputSet.iterItems(orderBy=['_micId', 'id'], direction='ASC'))

    for particle, row in zip(partIter, parFile):
        particle.setTransform(rowToAlignment(row, samplingRate))
        # We assume that each particle have ctfModel
        # in order to be processed in Frealign
        # JMRT: Since the CTF will be set, we can setup
        # an empty CTFModel object
        if not particle.hasCTF():
            particle.setCTF(CTFModel())
        rowToCtfModel(row, particle.getCTF())
        outputSet.append(particle)
    outputSet.setAlignment(ALIGN_PROJ)


def rowToCtfModel(ctfRow, ctfModel):
    """ Convert a row to Scipion CTF model.
    :param ctfRow: input row
    :param ctfModel: output model
    """
    defocusU = float(ctfRow.get('DF1'))
    defocusV = float(ctfRow.get('DF2'))
    defocusAngle = float(ctfRow.get('ANGAST'))
    ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)


def parseCtffindOutput(filename, avrot=False):
    """ Retrieve defocus U, V and angle from the
    output file of the ctffind execution.
    :param filename: input file to parse
    :param avrot: parse _avrot.txt file
    :return: an array with CTF values
    """
    if os.path.exists(filename):
        result = np.loadtxt(filename, dtype=float, comments='#')
        if avrot:
            # Reshape into a 3D array where each "block" of 6 lines forms one slice along the first axis
            reshaped_data = result.reshape(-1, 6, result.shape[1])
            # Apply the mask to the first line of each 6-line block
            # See details at https://github.com/3dem/relion-devel/commit/fceb6f687a09151dd60fd6e0ca46289798223bab
            mask = np.logical_and(reshaped_data[:, 0, :] >= 0.25, reshaped_data[:, 0, :] <= 0.28)
            # Compute the sum of the absolute values of the second line where the mask is True
            rotAvgArray = np.sum(np.abs(reshaped_data[:, 1, :]) * mask, axis=1)
            return rotAvgArray
        else:
            return result
    else:
        logger.error(f"Warning: Missing file: {filename}")

    return None


def setWrongDefocus(ctfModel):
    """ Set parameters if results parsing has failed.
    :param ctfModel: the model to be updated
    """
    ctfModel.setDefocusU(-999)
    ctfModel.setDefocusV(-1)
    ctfModel.setDefocusAngle(-999)
    
    
def readCtfModelStack(ctfModel, ctfArray, rotAvgArray=None, item=0):
    """ Set values for the ctfModel from an input list.
    :param ctfModel: output CTF model
    :param ctfArray: array with CTF values
    :param rotAvgArray: array with power spectrum values
    :param item: which row to use from ctfArray
    """
    if ctfArray is not None and ctfArray.ndim > 1:
        values = ctfArray[item]
    else:
        values = ctfArray

    tiltAxis, tiltAngle, thickness = None, None, None

    if ctfArray is None or np.isnan(values).any(axis=0) or values[1] < 0 or values[2] < 0:
        logger.debug(f"Invalid CTF values: {values}")
        setWrongDefocus(ctfModel)
        ctfFit, ctfResolution, ctfPhaseShift = -999, -999, 0
    else:
        defocusU, defocusV, defocusAngle = values[1], values[2], values[3]
        ctfPhaseShift, ctfFit, ctfResolution = values[4], values[5], values[6]
        ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)

        if len(values) == 10:
            tiltAxis, tiltAngle, thickness = values[7], values[8], values[9]

    ctfModel.setFitQuality(ctfFit)
    ctfModel.setResolution(ctfResolution)

    # Avoid creation of phaseShift
    ctfPhaseShiftDeg = np.rad2deg(ctfPhaseShift)
    if ctfPhaseShiftDeg != 0:
        ctfModel.setPhaseShift(ctfPhaseShiftDeg)

    if rotAvgArray is not None:
        if len(rotAvgArray) > 1:
            ctfModel._rlnIceRingDensity = Float(rotAvgArray[item])
        else:
            ctfModel._rlnIceRingDensity = Float(rotAvgArray[0])

    return ctfModel, tiltAxis, tiltAngle, thickness


def readCtfModel(ctfModel, ctfFn, rotAvgFn=None):
    """ Shortcut function to read CTFs for non-stacks
    and set ctfModel values.
    :param ctfModel: output CTF model
    :param ctfFn: input CTF file to parse
    :param rotAvgFn: extra input file with power spectra
    """
    result = parseCtffindOutput(ctfFn)

    rotAvgArray = None
    if rotAvgFn is not None:
        rotAvgArray = parseCtffindOutput(rotAvgFn, avrot=True)

    ctfModel, *_ = readCtfModelStack(ctfModel, ctfArray=result,
                                     rotAvgArray=rotAvgArray, item=0)

    return ctfModel


def readShiftsMovieAlignment(shiftFn):
    """ Parse movie alignment shifts.
    :param shiftFn: input file to parse
    :return: two lists with shift values
    """
    with open(shiftFn, 'r') as f:
        xshifts = []
        yshifts = []

        for line in f:
            line2 = line.strip()
            if line2.startswith('image #'):
                parts = line2.split()
                xshifts.append(float(parts[-2].rstrip(',')))
                yshifts.append(float(parts[-1]))
    return xshifts, yshifts


def readSetOfCoordinates(workDir, micSet, coordSet):
    """ Read coordinates from cisTEM .plt files.
    :param workDir: input folder with coord files
    :param micSet: input set of mics
    :param coordSet: output set of coords
    """
    for mic in micSet:
        micCoordFn = os.path.join(workDir, replaceBaseExt(mic.getFileName(), 'plt'))
        readCoordinates(mic, micCoordFn, coordSet)


def readCoordinates(mic, fn, coordsSet):
    """ Parse coords file and populate coordsSet.
    :param mic: input micrograph object
    :param fn: input file to parse
    :param coordsSet: output set of coords
    """
    if exists(fn):
        with open(fn, 'r') as f:
            for line in f:
                values = line.strip().split()
                # plt coords are in Imagic style
                x = float(values[1])
                y = float(mic.getYDim() - float(values[0]))
                coord = Coordinate()
                coord.setPosition(x, y)
                coord.setMicrograph(mic)
                coordsSet.append(coord)


def writeReferences(inputSet, outputFn):
    """ Write 2D references stack file from SetOfAverages or SetOfClasses2D.
    :param inputSet: the input SetOfParticles to be converted
    :param outputFn: where to write the output files.
    """
    ih = ImageHandler()

    def _convert(item, n):
        index = n + 1
        ih.convert(item, (index, outputFn))
        item.setLocation(index, outputFn)

    if isinstance(inputSet, SetOfAverages):
        for i, img in enumerate(inputSet):
            _convert(img, i)
    elif isinstance(inputSet, SetOfClasses2D):
        for i, rep in enumerate(inputSet.iterRepresentatives()):
            _convert(rep, i)
    else:
        raise TypeError('Invalid object type: %s' % type(inputSet))


def rowToAlignment(alignmentRow, samplingRate):
    """ Return an Transform object representing the Alignment
    from a given parFile row.
    :param alignmentRow: input row object
    :param samplingRate: input pixel size
    :return Transform object
    """
    angles = np.zeros(3)
    shifts = np.zeros(3)
    alignment = Transform()
    # PSI   THETA     PHI       SHX       SHY
    angles[0] = float(alignmentRow.get('PSI'))
    angles[1] = float(alignmentRow.get('THETA'))
    angles[2] = float(alignmentRow.get('PHI'))
    # shifts are converted from Angstroms to px
    shifts[0] = float(alignmentRow.get('SHX')) / samplingRate
    shifts[1] = float(alignmentRow.get('SHY')) / samplingRate

    M = matrixFromGeometry(shifts, angles)
    alignment.setMatrix(M)

    return alignment


def matrixFromGeometry(shifts, angles):
    """ Create the transformation matrix from given
    2D shifts in X and Y and the 3 euler angles.
    :param shifts: input list of shifts
    :param angles: input list of angles
    :return matrix
    """
    radAngles = -np.deg2rad(angles)

    M = transformations.euler_matrix(
        radAngles[0], radAngles[1], radAngles[2], 'szyz')
    M[:3, 3] = -shifts[:3]
    M = np.linalg.inv(M)

    return M


def geometryFromMatrix(matrix):
    """ Convert the transformation matrix to shifts and angles.
    :param matrix: input matrix
    :return: two lists, shifts and angles
    """
    matrix = np.linalg.inv(matrix)
    shifts = -transformations.translation_from_matrix(matrix)
    angles = -np.rad2deg(transformations.euler_from_matrix(matrix, axes='szyz'))

    return shifts, angles
