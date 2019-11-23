# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

import pyworkflow.em
import pyworkflow.utils as pwutils

from .constants import *


_logo = "cistem_logo.png"
_references = ['Grant2018']


class Plugin(pyworkflow.em.Plugin):
    _homeVar = CISTEM_HOME
    _pathVars = [CISTEM_HOME]
    _supportedVersions = [V1_0_0]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(CISTEM_HOME, 'cistem-1.0.0-beta')

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
        binary = os.path.join(cls.getHome(), program)
        if program == CTFFIND4_BIN:
            # if CTFFIND4_HOME is defined, use it
            path = os.environ.get(CTFFIND4_HOME, None) or cls.getHome()
            binary = os.path.join(path, program)

        return binary

    @classmethod
    def defineBinaries(cls, env):
        env.addPackage('cistem', version='1.0.0-beta',
                       url='"https://cistem.org/system/tdf/upload3/cistem-1.0.0'
                           '-beta-intel-linux.tar.gz?file=1&type=cistem_details'
                           '&id=37&force=0&s3fs=1"',
                       default=True)

        env.addPackage('ctffind4', version='4.1.13',
                       tar='ctffind4-4.1.13.tgz')


pyworkflow.em.Domain.registerPlugin(__name__)
