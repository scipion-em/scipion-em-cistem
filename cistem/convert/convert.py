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
from itertools import izip

from pyworkflow.em.data import Coordinate, SetOfClasses2D, SetOfAverages
from pyworkflow.em import ImageHandler
from pyworkflow.utils.path import replaceBaseExt, join, exists


def rowToCtfModel(ctfRow, ctfModel):
    defocusU = float(ctfRow.get('DF1'))
    defocusV = float(ctfRow.get('DF2'))
    defocusAngle = float(ctfRow.get('ANGAST'))
    ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)


def parseCtffind4Output(filename):
    """ Retrieve defocus U, V and angle from the
    output file of the ctffind4 execution.
    """
    result = None
    if os.path.exists(filename):
        f = open(filename)
        for line in f:
            if not line.startswith("#"):
                result = tuple(map(float, line.split()[1:]))
                # Stop reading. In ctffind4-4.0.15 output file has additional lines.
                break
        f.close()
    return result


def setWrongDefocus(ctfModel):
    ctfModel.setDefocusU(-999)
    ctfModel.setDefocusV(-1)
    ctfModel.setDefocusAngle(-999)
    
    
def readCtfModel(ctfModel, filename):
    result = parseCtffind4Output(filename)
    if result is None:
        setWrongDefocus(ctfModel)
        ctfFit, ctfResolution, ctfPhaseShift = -999, -999, -999
    else:
        defocusU, defocusV, defocusAngle, ctfPhaseShift, ctfFit, ctfResolution = result
        ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)
    ctfModel.setFitQuality(ctfFit)
    ctfModel.setResolution(ctfResolution)

    # Avoid creation of phaseShift
    ctfPhaseShiftDeg = np.rad2deg(ctfPhaseShift)
    if ctfPhaseShiftDeg != 0:
        ctfModel.setPhaseShift(ctfPhaseShiftDeg)


def readShiftsMovieAlignment(shiftFn):
    f = open(shiftFn, 'r')
    xshifts = []
    yshifts = []

    for line in f:
        l = line.strip()
        if l.startswith('image #'):
            parts = l.split()
            xshifts.append(float(parts[-2].rstrip(',')))
            yshifts.append(float(parts[-1]))
    f.close()
    return xshifts, yshifts


def writeShiftsMovieAlignment(movie, shiftsFn, s0, sN):
    movieAlignment = movie.getAlignment()
    shiftListX, shiftListY = movieAlignment.getShifts()

    # Generating metadata for global shifts
    a0, aN = movieAlignment.getRange()
    alFrame = a0

    if s0 < a0:
        diff = a0 - s0
        initShifts = "0.0000 " * diff
    else:
        initShifts = ""

    if sN > aN:
        diff = sN - aN
        finalShifts = "0.0000 " * diff
    else:
        finalShifts = ""

    shiftsX = ""
    shiftsY = ""
    for shiftX, shiftY in izip(shiftListX, shiftListY):
        if alFrame >= s0 and alFrame <= sN:
            shiftsX = shiftsX + "%0.4f " % shiftX
            shiftsY = shiftsY + "%0.4f " % shiftY
        alFrame += 1

    f = open(shiftsFn, 'w')
    shifts = (initShifts + shiftsX + " " + finalShifts + "\n"
              + initShifts + shiftsY + " " + finalShifts)
    f.write(shifts)
    f.close()


def readSetOfCoordinates(workDir, micSet, coordSet):
    """ Read from cisTEM .plt files. """
    for mic in micSet:
        micCoordFn = join(workDir, replaceBaseExt(mic.getFileName(), 'plt'))
        readCoordinates(mic, micCoordFn, coordSet)


def readCoordinates(mic, fn, coordsSet):
    if exists(fn):
        with open(fn, 'r') as f:
            for line in f:
                values = line.strip().split()
                x = float(values[1])
                y = float(mic.getYDim() - float(values[0]))
                coord = Coordinate()
                coord.setPosition(x, y)
                coord.setMicrograph(mic)
                coordsSet.append(coord)
        f.close()


def writeReferences(inputSet, outputFn):
    """
    Write references star and stack files from SetOfAverages or SetOfClasses2D.
    Params:
        inputSet: the input SetOfParticles to be converted
        outputFn: where to write the output files.
    """
    ih = ImageHandler()

    def _convert(item, i):
        index = i + 1
        ih.convert(item, (index, outputFn))
        item.setLocation(index, outputFn)

    if isinstance(inputSet, SetOfAverages):
        for i, img in enumerate(inputSet):
            _convert(img, i)
    elif isinstance(inputSet, SetOfClasses2D):
        for i, rep in enumerate(inputSet.iterRepresentatives()):
            _convert(rep, i)
    else:
        raise Exception('Invalid object type: %s' % type(inputSet))
