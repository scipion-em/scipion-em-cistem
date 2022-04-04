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

import os
from matplotlib.figure import Figure

from pyworkflow.utils import removeExt
from pwem.viewers import EmPlotter
from pwem.emlib.image import ImageHandler
from tomo.viewers.viewers_data import CtfEstimationTomoViewer
from cistem.protocols import CistemProtTsCtffind
from cistem.viewers.viewers import getPlotSubtitle, _getValues


class CtfEstimationTomoViewerCistem(CtfEstimationTomoViewer):
    """ This class implements a view using Tkinter CtfEstimationListDialog
    and the CtfEstimationTreeProvider.
    """
    _targets = [CistemProtTsCtffind]

    def plot1D(self, ctfSet, ctfId):
        ctfModel = ctfSet[ctfId]
        psdFn = ctfModel.getPsdFile()
        fn = os.path.join(removeExt(psdFn) + '_avrot.txt')

        xplotter = EmPlotter(windowTitle='CTFFind results')
        plot_title = '%s # %d\n' % (ctfSet.getTsId(), ctfId) + getPlotSubtitle(ctfModel)
        a = xplotter.createSubPlot(plot_title, 'Spacial frequency (1/A)',
                                   'Amplitude (or cross-correlation)')
        legendName = ['Amplitude spectrum',
                      'CTF Fit',
                      'Quality of fit']
        res = _getValues(fn)
        for y in ['amp', 'fit', 'quality']:
            a.plot(res['freq'], res[y])
        xplotter.showLegend(legendName, loc='upper right')
        a.set_ylim([-0.1, 1.1])
        a.grid(True)

        return xplotter

    def plot2D(self, ctfSet, ctfId):
        ctfModel = ctfSet[ctfId]
        psdFn = ctfModel.getPsdFile()
        if not os.path.exists(psdFn):
            return None
        img = ImageHandler().read(psdFn)
        fig = Figure(figsize=(7, 7), dpi=100)
        psdPlot = fig.add_subplot(111)
        psdPlot.get_xaxis().set_visible(False)
        psdPlot.get_yaxis().set_visible(False)
        psdPlot.set_title('%s # %d\n' % (ctfSet.getTsId(), ctfId) + getPlotSubtitle(ctfModel))
        psdPlot.imshow(img.getData(), cmap='gray')

        return fig
