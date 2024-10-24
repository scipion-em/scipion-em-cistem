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
import re

import pwem
import pyworkflow.utils as pwutils

from .constants import *


__version__ = '3.10'
_logo = "cistem_logo.png"
_references = ['Grant2018']


class Plugin(pwem.Plugin):
    _homeVar = CISTEM_HOME
    _pathVars = [CISTEM_HOME, CTFFIND_HOME]
    _supportedVersions = [V1_0_0]
    _url = "https://github.com/scipion-em/scipion-em-cistem"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(CISTEM_HOME, 'cistem-1.0.0-beta')
        cls._defineEmVar(CTFFIND_HOME, 'ctffind-5.0.2')

    @classmethod
    def getActiveVersion(cls, *args):
        """ Return the env name that is currently active. """
        ctffind = cls.getVar(CTFFIND_HOME)
        return re.search(r"ctffind-([0-9a-zA-Z.]+)$", ctffind).group(1)

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
        if program == CTFFIND_BIN:
            # if CTFFIND_HOME is found, use it
            path = cls.getVar(CTFFIND_HOME)
            if os.path.exists(path):
                binary = os.path.join(path, 'bin', program)
            else:
                binary = os.path.join(cls.getHome(), program)
        else:
            binary = os.path.join(cls.getHome(), program)

        return binary

    @classmethod
    def defineBinaries(cls, env):
        env.addPackage('cistem', version='1.0.0-beta',
                       url="https://grigoriefflab.umassmed.edu/sites/default/files/cistem-1.0.0-beta-intel-linux.tar.gz",
                       default=True)
        env.addPackage('ctffind4', version='4.1.14',
                       tar='ctffind4-4.1.14.tgz')
        env.addPackage('ctffind', version='5.0',
                       tar='ctffind-5.0.tgz')
        env.addPackage('ctffind', version='5.0.2',
                       tar='ctffind-5.0.2.tgz',
                       default=True)
