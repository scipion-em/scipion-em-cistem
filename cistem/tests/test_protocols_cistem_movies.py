# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *             Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from cistem.protocols import CistemProtUnblur



class TestMoviesBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('movies')
        cls.movies = cls.dataset.getFile('ribo/*.mrcs')
    
    @classmethod
    def runImportMovie(cls, pattern, samplingRate, voltage,
                       sphericalAberration):
        """ Run an Import movies protocol. """
        kwargs = {
            'filesPath': pattern,
            'voltage': voltage,
            'sphericalAberration': sphericalAberration,
            'amplitudeContrast': 0.1,
            'dosePerFrame': 1.3,
            'samplingRate': samplingRate
        }

        cls.protImport = cls.newProtocol(ProtImportMovies, **kwargs)
        cls.proj.launchProtocol(cls.protImport, wait=True)

        return cls.protImport

    @classmethod
    def runImportMovies(cls):
        """ Run an Import movie protocol. """
        return cls.runImportMovie(cls.movies,
                                  samplingRate=3.54, voltage=300,
                                  sphericalAberration=2.0)


class TestUnblur(TestMoviesBase):
    def test_movies(self):
        protImport = self.runImportMovies()
        outputMovies = getattr(protImport, 'outputMovies', None)
        self.assertIsNotNone(outputMovies)

        inputMovies = protImport.outputMovies
        prot = self.newProtocol(CistemProtUnblur,
                                alignFrame0=2,
                                alignFrameN=6)
        prot.inputMovies.set(inputMovies)
        self.launchProtocol(prot)

        outputMics = getattr(prot, 'outputMicrographsDoseWeighted', None)
        self.assertIsNotNone(outputMics)
        self.assertEqual(protImport.outputMovies.getSize(),
                         outputMics.getSize())

        for mic in outputMics:
            micFn = mic.getFileName()
            self.assertTrue(os.path.exists(self.proj.getPath(micFn)))
