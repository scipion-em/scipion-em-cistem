# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow.utils as pwutils 
from pwem.objects import CTFModel, SetOfParticles

from .convert import readCtfModel, readSetOfParticles


class GrigorieffLabImportCTF:
    """ Import CTF estimated with CTFFIND4. """
    def __init__(self, protocol):
        self.protocol = protocol
        self.copyOrLink = self.protocol.getCopyOrLink()

    def importCTF(self, mic, fileName):
        """ Create a CTF model and populate its values.
        :param mic: input micrograph object
        :param fileName: input file to be parsed
        :return: CTFModel object
        """
        ctf = CTFModel()
        ctf.setMicrograph(mic)
        readCtfModel(ctf, fileName)
        
        # Try to find the given PSD file associated with the cttfind log file
        # we handle special cases of .ctf extension and _ctffind4 prefix for Relion runs
        fnBase = pwutils.removeExt(fileName)
        for suffix in ['_psd.mrc', '.mrc', '.ctf']:
            psdPrefixes = [fnBase,
                           fnBase.replace('_ctffind4', '')]
            for prefix in psdPrefixes:
                psdFile = prefix + suffix
                if os.path.exists(psdFile):
                    if psdFile.endswith('.ctf'):
                        psdFile += ':mrc'
                    ctf.setPsdFile(psdFile)
                    return ctf

        return ctf


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
