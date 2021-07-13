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

import contextlib

from .protocol_ctffind import CistemProtCTFFind
from .protocol_unblur import CistemProtUnblur
from .protocol_picking import CistemProtFindParticles
from .protocol_refine2d import CistemProtRefine2D

# This is a prototype and will likely be move to pyworkflow soon
# Once in pyworkflow, this method can be removed and imported from there
@contextlib.contextmanager
def weakImport(package):
    """
     This method can be use to tolerate imports that may fail, e.g imports
    :param package: name of the package that is expected to fail
    """
    try:
        yield
    except ImportError as e:
        if "'%s'" % package not in str(e):
            raise e


with weakImport('tomo'):
    from .protocol_ts_ctffind import CistemProtTsCtffind

