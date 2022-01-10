# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] Department of Anatomy and Cell Biology, McGill University
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

from pyworkflow.protocol.params import LabelParam
from pyworkflow.utils import removeExt, cleanPath
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer
from pyworkflow.gui.project import ProjectWindow
from pwem.viewers import CtfView, EmPlotter, MicrographsView, EmProtocolViewer
import pwem.viewers.showj as showj
from pwem.objects import SetOfMovies


from cistem.protocols import CistemProtCTFFind, CistemProtUnblur


def createCtfPlot(ctfSet, ctfId):
    """ Create EmPlotter instance. """
    ctfModel = ctfSet[ctfId]
    psdFn = ctfModel.getPsdFile()
    fn = removeExt(psdFn) + "_avrot.txt"
    xplotter = EmPlotter(windowTitle='CTFFind results')
    plot_title = getPlotSubtitle(ctfModel)
    a = xplotter.createSubPlot(plot_title, 'Spacial frequency (1/A)',
                               'Amplitude (or cross-correlation)')
    legendName = ['Amplitude spectrum',
                  'CTF Fit',
                  'Quality of fit']
    _plotCurves(a, fn)
    xplotter.showLegend(legendName, loc='upper right')
    a.set_ylim([-0.1, 1.1])
    a.grid(True)
    xplotter.show()


def getPlotSubtitle(ctf):
    """ Create plot subtitle using CTF values. """
    ang = u"\u212B"
    deg = u"\u00b0"
    def1, def2, angle = ctf.getDefocus()
    phSh = ctf.getPhaseShift()
    score = ctf.getFitQuality()
    res = ctf.getResolution()

    title = "Def1: %d %s | Def2: %d %s | Angle: %0.1f%s | " % (
        def1, ang, def2, ang, angle, deg)

    if phSh is not None:
        title += "Phase shift: %0.2f %s | " % (phSh, deg)

    title += "Fit: %0.1f %s | Score: %0.3f" % (res, ang, score)

    return title


def _plotCurves(a, fn):
    """ Actually plot the curves. """
    res = _getValues(fn)
    for y in ['amp', 'fit', 'quality']:
        a.plot(res['freq'], res[y])


def _getValues(fn):
    """ Parse input file and return a dict with results. """
    res = dict()
    with open(fn) as f:
        i = 0
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                if i == 0:
                    res['freq'] = [float(x) for x in line.split()]
                elif i == 2:
                    res['amp'] = [float(x) for x in line.split()]
                elif i == 3:
                    res['fit'] = [float(x) for x in line.split()]
                elif i == 4:
                    res['quality'] = [float(x) for x in line.split()]
                    break
                i += 1
    return res


OBJCMD_CTFFIND4 = "CTFFind plot results"

ProjectWindow.registerObjectCommand(OBJCMD_CTFFIND4, createCtfPlot)


class CtffindViewer(Viewer):
    """ Specific way to visualize SetOfCtf. """
    _environments = [DESKTOP_TKINTER]
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


class ProtUnblurViewer(EmProtocolViewer):
    _targets = [CistemProtUnblur]
    _environments = [DESKTOP_TKINTER]

    _label = 'viewer unblur'

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        if self.hasMics():
            form.addParam('doShowMics', LabelParam,
                          label="Show aligned micrographs?", default=True,
                          help="Show the output aligned micrographs.")

        if self.hasDWMics():
            form.addParam('doShowMicsDW', LabelParam,
                          label="Show aligned DOSE-WEIGHTED micrographs?",
                          default=True,
                          help="Show the output aligned dose-weighted "
                               "micrographs.")
        form.addParam('doShowMovies', LabelParam,
                      label="Show output movies?", default=True,
                      help="Show the output movies with alignment "
                           "information.")
        form.addParam('doShowFailedMovies', LabelParam,
                      label="Show FAILED movies?", default=True,
                      help="Create a set of failed movies "
                           "and display it.")

    def _getVisualizeDict(self):
        self._errors = []

        visualizeDict = {
                         'doShowMovies': self._viewParam,
                         'doShowFailedMovies': self._viewParam
                         }
        if self.hasMics():
            visualizeDict.update({'doShowMics': self._viewParam})

        if self.hasDWMics():
            visualizeDict.update({'doShowMicsDW': self._viewParam})

        return visualizeDict

    def hasMics(self):
        return hasattr(self.protocol, 'outputMicrographs')

    def hasDWMics(self):
        return hasattr(self.protocol, 'outputMicrographsDoseWeighted')

    def _viewParam(self, param=None):
        labelsDef = 'enabled id _filename _samplingRate '
        labelsDef += '_acquisition._dosePerFrame _acquisition._doseInitial '
        viewParamsDef = {showj.MODE: showj.MODE_MD,
                         showj.ORDER: labelsDef,
                         showj.VISIBLE: labelsDef,
                         showj.RENDER: None
                         }
        if param == 'doShowMics':
            return [MicrographsView(self.getProject(),
                                    self.protocol.outputMicrographs)]
        elif param == 'doShowMicsDW':
            return [MicrographsView(self.getProject(),
                                    self.protocol.outputMicrographsDoseWeighted)]

        elif param == 'doShowMovies':
            if getattr(self.protocol, 'outputMovies', None) is not None:
                output = self.protocol.outputMovies
                return [self.objectView(output, viewParams=viewParamsDef)]
            else:
                return [self.errorMessage('No output movies found!',
                                          title="Visualization error")]

        elif param == 'doShowFailedMovies':
            self.failedList = self.protocol._readFailedList()
            if not self.failedList:
                return [self.errorMessage('No failed movies found!',
                                          title="Visualization error")]
            else:
                sqliteFn = self.protocol._getPath('movies_failed.sqlite')
                self.createFailedMoviesSqlite(sqliteFn)
                return [self.objectView(sqliteFn, viewParams=viewParamsDef)]

    def createFailedMoviesSqlite(self, path):
        inputMovies = self.protocol.inputMovies.get()
        cleanPath(path)
        movieSet = SetOfMovies(filename=path)
        movieSet.copyInfo(inputMovies)
        movieSet.copyItems(inputMovies,
                           updateItemCallback=self._findFailedMovies)
        movieSet.write()
        movieSet.close()

        return movieSet

    def _findFailedMovies(self, item, row):
        if item.getObjId() not in self.failedList:
            setattr(item, "_appendItem", False)
