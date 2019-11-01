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

from pyworkflow.utils import removeExt
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewers import CtfView, EmPlotter, MicrographsView
import pyworkflow.em.viewers.showj as showj
from pyworkflow.gui.project import ProjectWindow

from cistem.protocols import CistemProtCTFFind, CistemProtUnblur


def createCtfPlot(ctfSet, ctfId):
    ctfModel = ctfSet[ctfId]
    psdFn = ctfModel.getPsdFile()
    fn = removeExt(psdFn) + "_avrot.txt"
    gridsize = [1, 1]
    xplotter = EmPlotter(x=gridsize[0], y=gridsize[1],
                         windowTitle='CTF Fitting')
    plot_title = "CTF Fitting"
    a = xplotter.createSubPlot(plot_title, 'pixels^-1', 'CTF',
                               yformat=False)
    
    legendName = ['rotational avg. No Astg',
                  'rotational avg.',
                  'CTF Fit',
                  'Cross Correlation',
                  '2sigma cross correlation of noise']
    for i in range(1, 6):
        _plotCurve(a, i, fn)
    xplotter.showLegend(legendName)
    a.grid(True)
    xplotter.show()


OBJCMD_CTFFIND4 = "Display Ctf Fitting"

ProjectWindow.registerObjectCommand(OBJCMD_CTFFIND4, createCtfPlot)


class CtffindViewer(Viewer):
    """ Specific way to visualize SetOfCtf. """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [CistemProtCTFFind]

    def _visualize(self, prot, **kwargs):
        outputCTF = getattr(prot, 'outputCTF', None)

        if outputCTF is not None:
            ctfView = CtfView(self._project, outputCTF)
            viewParams = ctfView.getViewParams()
            viewParams[showj.OBJCMDS] = "'%s'" % OBJCMD_CTFFIND4
            return [ctfView]
        else:
            return [self.infoMessage("The output SetOfCTFs has not been "
                                     "produced", "Missing output")]


def _plotCurve(a, i, fn):
    freqs = _getValues(fn, 0)
    curv = _getValues(fn, i)
    a.plot(freqs, curv)

def _getValues(fn, row):
    f = open(fn)
    values = []
    i = 0
    for line in f:
        if not line.startswith("#"):
            if i == row:
                values = line.split()
                break
            i += 1
    f.close()
    return values


class ProtUnblurViewer(Viewer):
    _targets = [CistemProtUnblur]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    _label = 'viewer unblur'

    def _visualize(self, obj, **kwargs):
        views = []

        labelsDef = 'enabled id _filename'
        viewParamsDef = {showj.MODE: showj.MODE_MD,
                         showj.ORDER: labelsDef,
                         showj.VISIBLE: labelsDef,
                         showj.RENDER: None
                         }

        outputLabels = ['outputMicrographs', 'outputMicrographsDoseWeighted',
                        'outputMovies']

        if not any(obj.hasAttribute(l) for l in outputLabels):
            return [self.infoMessage("Output (micrographs or movies) have "
                                     "not been produced yet.")]

        # Display only the first available output, showing all of them
        # can be confusing and not useful.
        # The user can still double-click in the specific output
        for l in outputLabels:
            if obj.hasAttribute(l):
                output = getattr(obj, l)
                if 'Micrographs' in l:
                    return [MicrographsView(self.getProject(), output)]
                else:  # Movies case
                    return [self.objectView(output, viewParams=viewParamsDef)]

        return views
