# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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

import pyworkflow.em
import pyworkflow.utils as pwutils

from .constants import CISTEM_HOME, V1_0_0


_logo = "cistem_logo.png"
_references = ['Grant2018']



class Plugin(pyworkflow.em.Plugin):
    _homeVar = CISTEM_HOME
    _pathVars = [CISTEM_HOME]
    _supportedVersions = [V1_0_0]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(CISTEM_HOME, 'cistem-1.0.0')

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch cisTEM. """
        environ = pwutils.Environ(os.environ)
        environ.update({'PATH': cls.getHome()},
                       position=pwutils.Environ.BEGIN)

        return environ

    @classmethod
    def getProgram(cls, program):
        """ Return the program binary that will be used. """
        program = os.path.join(cls.getHome(), program)
        return program

    @classmethod
    def defineBinaries(cls, env):
        env.addPackage('cistem', version='1.0.0',
                       tar='cistem-1.0.0-beta-intel-linux.tar.gz',
                       default=True)


pyworkflow.em.Domain.registerPlugin(__name__)
