# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] SciLifeLab, Stockholm University
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

import os.path

import pyworkflow.utils as pwutils
from pwem.objects import CTFModel, SetOfParticles

from .convert import (readCtfModel, readSetOfParticles,
                      readCtfModelStack, parseCtffindOutput)

with pwutils.weakImport('tomo'):
    from tomo.objects import CTFTomo


class GrigorieffLabImportCTF:
    """ Import CTF estimated with CTFFIND. """
    def __init__(self, protocol):
        self.protocol = protocol
        self.copyOrLink = self.protocol.getCopyOrLink()

    def importCTF(self, mic, fileName):
        """ Create a CTF model and populate its values.
        :param mic: input micrograph object
        :param fileName: input file to be parsed
        :return: CTFModel object
        """
        fnBase = pwutils.removeExt(fileName)
        avrotFn = fileName + "_avrot.txt"
        avrotFn = None if not os.path.exists(avrotFn) else avrotFn

        ctf = CTFModel()
        ctf.setMicrograph(mic)
        ctf = readCtfModel(ctf, fileName, avrotFn)

        psdFile = self._findPsdFile(fnBase)
        ctf.setPsdFile(psdFile)

        return ctf

    def parseTSDefocusFile(self, ts, fileName, output):
        """ Parse tilt-series ctf estimation file.
        :param ts: input tilt-series
        :param fileName: input file to be parsed
        :param output: output CTFTomoSeries
        """
        tsId = ts.getTsId()
        fnBase = os.path.join(os.path.dirname(fileName), tsId)
        outputPsd = self._findPsdFile(fnBase)
        ctfResult = parseCtffindOutput(fileName)
        avrotFn = pwutils.removeExt(fileName) + "_avrot.txt"
        if os.path.exists(avrotFn):
            avrotResult = parseCtffindOutput(avrotFn, avrot=True)
        else:
            avrotResult = None

        ctf = CTFModel()
        counter = 0

        for i, ti in enumerate(ts):
            if ti.isEnabled():
                self.getCtfTi(ctf, ctfResult, avrotResult, counter, outputPsd)
                counter += 1
            else:
                ctf.setWrongDefocus()

            newCtfTomo = CTFTomo.ctfModelToCtfTomo(ctf)
            if hasattr(ctf, "_rlnIceRingDensity"):
                newCtfTomo.copyAttributes(ctf, '_rlnIceRingDensity')
            if not ti.isEnabled():
                newCtfTomo.setEnabled(False)
            newCtfTomo.setIndex(i + 1)
            newCtfTomo.setAcquisitionOrder(ti.getAcquisitionOrder())
            output.append(newCtfTomo)

    @staticmethod
    def getCtfTi(ctf, ctfArray, rotAvgArray, tiIndex, psdStack=None):
        """ Parse the CTF object estimated for this Tilt-Image. """
        ctf, tiltAxis, tiltAngle, thickness = readCtfModelStack(ctf, ctfArray, rotAvgArray, item=tiIndex)
        if psdStack is not None:
            ctf.setPsdFile(f"{tiIndex + 1}@" + psdStack)

    @staticmethod
    def _findPsdFile(fnBase):
        """ Try to find the given PSD file associated with the cttfind log file
        We handle special cases of .ctf extension and _ctffind4 prefix for Relion runs
        """
        for suffix in ['_psd.mrc', '.mrc', '_ctf.mrcs',
                       '.mrcs', '.ctf']:
            psdPrefixes = [fnBase,
                           fnBase.replace('_ctffind4', '')]
            for prefix in psdPrefixes:
                psdFile = prefix + suffix
                if os.path.exists(psdFile):
                    if psdFile.endswith('.ctf'):
                        psdFile += ':mrc'
                    return psdFile
        return None


class GrigorieffLabImportParticles:
    """ Import particles from a Frealign refinement.
    :param parFile: the filename of the parameter file with the alignment
    :param stackFile: single stack file with the images
    """
    def __init__(self, protocol, parFile, stackFile):
        self.protocol = protocol
        self.copyOrLink = self.protocol.getCopyOrLink()
        self.parFile = parFile
        self.stackFile = stackFile

    def _setupSet(self, partSet):
        self.protocol.setSamplingRate(partSet)
        partSet.setIsPhaseFlipped(self.protocol.haveDataBeenPhaseFlipped.get())
        self.protocol.fillAcquisition(partSet.getAcquisition())

    def importParticles(self):
        partSet = self.protocol._createSetOfParticles()
        partSet.setObjComment('Particles imported from Frealign parfile:\n%s' % self.parFile)

        # Create a local link to the input stack file
        localStack = self.protocol._getExtraPath(os.path.basename(self.stackFile))
        pwutils.createLink(self.stackFile, localStack)
        # Create a temporary set only with location
        tmpSet = SetOfParticles(filename=':memory:')
        tmpSet.readStack(localStack)
        self._setupSet(tmpSet)

        # Update both samplingRate and acquisition with parameters
        # selected in the protocol form
        self._setupSet(partSet)
        # Now read the alignment parameters from par file
        readSetOfParticles(tmpSet, partSet, self.parFile)
        partSet.setHasCTF(True)
        # Register the output set of particles
        self.protocol._defineOutputs(outputParticles=partSet)

    def validateParticles(self):
        """ Overwrite the base class. """
        errors = []
        return errors
