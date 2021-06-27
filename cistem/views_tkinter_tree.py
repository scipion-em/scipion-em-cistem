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

from tkinter import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from pyworkflow.gui import *
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.gui.dialog import ListDialog, showInfo

import tomo.objects


class CTFSerieStates:
    UNCHECKED = 'unchecked'
    CHECKED = 'checked'
    ODD = 'odd'
    EVEN = 'even'
    FAILED = 'Failed'
    OK = 'Ok'


class CtfEstimationTreeProvider(TreeProvider, ttk.Treeview):
    """ Model class that will retrieve the information from TiltSeries and
    prepare the columns/rows models required by the TreeDialog GUI.
    """
    COL_CTF_SERIE = 'Tilt Series'
    COL_TILT_ANG = 'Tilt Angle'
    CRITERIA_1 = 'Status'
    COL_CTF_EST_IX = 'Index'
    COL_CTF_EST_DEFOCUS_U = 'Defocus (A)'
    COL_CTF_EST_AST = 'Astigmatism (A)'
    COL_CTF_EST_RES = 'Resolution (A)'
    COL_CTF_EST_FIT = 'CC value'
    COL_CTF_EST_PHASE = 'Phase shift (deg)'

    ORDER_DICT = {COL_CTF_EST_DEFOCUS_U: '_defocusU',
                  COL_CTF_EST_AST: '_defocusRatio',
                  COL_CTF_EST_RES: '_resolution',
                  COL_CTF_EST_PHASE: '_phaseShift',
                  COL_CTF_EST_FIT: '_fitQuality'}

    def __init__(self, master, protocol, outputSetOfCTFTomoSeries, **kw):
        ttk.Treeview.__init__(self, master, **kw)
        self.protocol = protocol
        self.ctfSeries = outputSetOfCTFTomoSeries
        TreeProvider.__init__(self, sortingColumnName=self.COL_CTF_SERIE)
        self.selectedDict = {}
        self.mapper = protocol.mapper
        self.maxNum = 200
        self._checkedItems = 0

    def getObjects(self):
        # Retrieve all objects of type className
        objects = []

        orderBy = self.ORDER_DICT.get(self.getSortingColumnName(), 'id')
        direction = 'ASC' if self.isSortingAscending() else 'DESC'

        for ctfSerie in self.ctfSeries:
            ctfEstObj = ctfSerie.clone()
            ctfEstObj._allowsSelection = True
            ctfEstObj._parentObject = None
            objects.append(ctfEstObj)
            for item in ctfSerie.iterItems(orderBy=orderBy, direction=direction):
                ctfEstItem = item.clone()
                ctfEstItem._allowsSelection = False
                ctfEstItem._parentObject = ctfEstObj
                objects.append(ctfEstItem)

        return objects

    def getCTFSeries(self):
        return self.ctfSeries

    def _sortObjects(self, objects):
        pass

    def objectKey(self, pobj):
        pass

    def showPhaseShiftCol(self):
        ctfSerie = self.ctfSeries.getFirstItem()
        return ctfSerie.getFirstItem().hasPhaseShift()

    def getColumns(self):
        cols = [
            (self.COL_CTF_SERIE, 100),
            (self.COL_TILT_ANG, 100),
            (self.CRITERIA_1, 60),
            (self.COL_CTF_EST_DEFOCUS_U, 100),
            (self.COL_CTF_EST_AST, 150),
        ]
        if self.showPhaseShiftCol():
            cols.extend([
                (self.COL_CTF_EST_PHASE, 150),
                (self.COL_CTF_EST_RES, 100),
                (self.COL_CTF_EST_FIT, 100)
            ])
        else:
            cols.extend([
                (self.COL_CTF_EST_RES, 100),
                (self.COL_CTF_EST_FIT, 100)
            ])

        return cols

    def isSelected(self, obj):
        """ Check if an object is selected or not. """
        return False

    @staticmethod
    def _getParentObject(pobj, default=None):
        return getattr(pobj, '_parentObject', default)

    def getObjectInfo(self, obj):
        if isinstance(obj, tomo.objects.CTFTomoSeries):
            key = obj.getTsId()
            text = obj.getTsId()
            values = ['', CTFSerieStates.OK if obj.getIsDefocusUDeviationInRange()
                      else CTFSerieStates.FAILED]
            opened = False
            selected = obj.isEnabled()
        else:  # CTFTomo
            key = "%s.%s" % (obj._parentObject.getTsId(), str(obj.getObjId()))
            text = obj.getIndex()
            ts = obj._parentObject.getTiltSeries()
            tiltAngle = ts[int(text)].getTiltAngle()
            ast = obj.getDefocusU() - obj.getDefocusV()
            phSh = obj.getPhaseShift() if obj.hasPhaseShift() else 0

            values = [str("%0.2f" % tiltAngle),
                      CTFSerieStates.OK if obj.getIsDefocusUDeviationInRange()
                      else CTFSerieStates.FAILED,
                      str("%d" % obj.getDefocusU()),
                      str("%d" % ast)]

            if self.showPhaseShiftCol():
                values.extend([
                    str("%0.2f" % phSh),
                    str("%0.1f" % obj.getResolution()),
                    str("%0.2f" % obj.getFitQuality())
                ])
            else:
                values.extend([
                    str("%0.1f" % obj.getResolution()),
                    str("%0.2f" % obj.getFitQuality())
                ])

            opened = False
            selected = False

        item = {
            'key': key, 'text': text,
            'values': tuple(values),
            'open': opened,
            'selected': selected,
            'parent': obj._parentObject
        }
        if isinstance(obj, tomo.objects.CTFTomoSeries):
            tags = CTFSerieStates.UNCHECKED
            if not (obj.getIsDefocusUDeviationInRange() and obj.getIsDefocusVDeviationInRange()):
                obj.setEnabled(True)
                tags = CTFSerieStates.CHECKED
                self._checkedItems += 1

            if obj.getObjId() % 2 == 0:
                item['tags'] = (tags, CTFSerieStates.ODD)
            else:
                item['tags'] = (tags,  CTFSerieStates.EVEN)
        else:
            if obj.getObjId() % 2 == 0:
                item['tags'] = (CTFSerieStates.ODD,)
            else:
                item['tags'] = (CTFSerieStates.EVEN,)
        return item


class CTFEstimationTree(BoundTree):
    def __init__(self, master, provider,  **opts):
        BoundTree.__init__(self, master, provider, **opts)
        self.selectedItem = None
        self._checkedItems = provider._checkedItems

    def check_item(self, item):
        """ check the box of item and change the state of the boxes of item's
            ancestors accordingly """
        tags = CTFSerieStates.EVEN
        if CTFSerieStates.ODD in self.item(item, 'tags'):
            tags = CTFSerieStates.ODD

        if CTFSerieStates.UNCHECKED in self.item(item, 'tags'):
            self.item(item, tags=(CTFSerieStates.CHECKED, tags,))
            self._checkedItems += 1
            self.getSelectedObj().setEnabled(False)
            self.item(item)['selected'] = True
        else:
            self.item(item, tags=(CTFSerieStates.UNCHECKED, tags,))
            self.getSelectedObj().setEnabled(True)
            self._checkedItems -= 1
            self.item(item)['selected'] = False

    def _onClick(self, event=None):
        self._unpostMenu()
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        self.selectedItem = self.identify_row(y)
        self.focus(self.selectedItem)
        if "image" in elem:  # click on the checkbox
            self.check_item(self.selectedItem)

    def getSelectedItem(self):
        return self.selectedItem

    def getSelectedObj(self):
        obj = None
        if self.selectedItem:
            obj = self._objDict[self.getFirst()]
        return obj


class CtfEstimationListDialog(ListDialog):
    def __init__(self, parent, title, provider, protocol, inputTS, **kwargs):
        self._project = protocol.getProject()
        self._protocol = protocol
        self._inputSetOfTiltSeries = inputTS
        self._checkedItems = provider._checkedItems
        # the func below should be implemented in viewers
        self._show1DPLot = kwargs.pop('plot1Dfunc')
        self._show2DPLot = kwargs.pop('plot2Dfunc')
        ListDialog.__init__(self, parent, title, provider, allowSelect=False,
                            cancelButton=True, **kwargs)

    def body(self, bodyFrame):
        bodyFrame.config()
        self._col = 1
        self._fillCTFEstimationGUI(bodyFrame)

    def _addButton(self, frame, text, image, command, sticky='news', state=tk.NORMAL):
        btn = tk.Label(frame, text=text, image=self.getImage(image),
                       compound=tk.LEFT, cursor='hand2', state=state)
        btn.bind('<Button-1>', command)
        btn.grid(row=0, column=self._col, sticky=sticky,
                 padx=(0, 5), pady=5)
        self._col += 1
        return btn

    def _fillCTFEstimationGUI(self, bodyFrame):
        # Create a top panel to put the filter box and bottoms
        topPanel = tk.Frame(bodyFrame)
        topPanel.grid(row=0, column=0, padx=0, pady=0, sticky='news')
        self._createTopPanel(topPanel)

        # Create a bottom panel to put the tree and the plotter
        bottomPanel = tk.Frame(bodyFrame)
        bottomPanel.grid(row=1, column=0, padx=0, pady=0, sticky='news')
        self._createBottomPanel(bottomPanel)

    def _createTopPanel(self, topPanel):
        self._createFilterBox(topPanel)

        topRigthPanel = tk.Frame(topPanel)
        topRigthPanel.grid(row=0, column=1, padx=0, pady=0, sticky='news')
        self._createRecalculateBottom(topRigthPanel)
        self._createShowFit(topRigthPanel)
        self._createViewerHelp(topRigthPanel)

    def _createRecalculateBottom(self, topRigthPanel):
        state = tk.NORMAL
        if self._checkedItems or self._checkedItems == len(self.provider.getCTFSeries()):
            state = tk.DISABLED
        self.generateSubsetButton = self._addButton(topRigthPanel,
                                                    'Generate subsets',
                                                    pwutils.Icon.PROCESSING,
                                                    self._actionCreateSets,
                                                    sticky='ne',
                                                    state=state)

    def _createShowFit(self, topRigthPanel):
        self._addButton(topRigthPanel, '1D fit',
                        pwutils.Icon.ACTION_RESULTS, self._show1DFit, sticky='ne')

    def _createViewerHelp(self, topRigthPanel):
        self._addButton(topRigthPanel, pwutils.Message.LABEL_HELP,
                        pwutils.Icon.ACTION_HELP, self._showHelp, sticky='ne')

    def _show1DFit(self, event=None):
        itemSelected = self.tree.getSelectedItem()
        obj = self.tree.getSelectedObj()
        if self.tree.parent(itemSelected):  # child item
            if obj is not None:
                for ctfSerie in self.provider.getCTFSeries():
                    if ctfSerie.getTsId() in itemSelected:
                        # TODO: is ctfSerie ordered by id?
                        ctfId = int(itemSelected.split('.')[-1])
                        plot = self._show1DPLot(ctfSerie, ctfId)
                        plot.show()
                        break

    def _actionCreateSets(self, event=None):
        if self.generateSubsetButton['state'] == tk.NORMAL:
            protocol = self.provider.protocol
            ctfSeries = self.provider.getCTFSeries()
            suffix = self._getSuffix(protocol)
            goodCTFName = 'goodCtf%s' % suffix
            badCTFName = 'badCtf%s' % suffix

            outputSetOfgoodCTFTomoSeries = ctfSeries.createCopy(protocol._getPath(),
                                                                prefix=goodCTFName,
                                                                copyInfo=True)
            outputSetOfbadCTFTomoSeries = ctfSeries.createCopy(protocol._getPath(),
                                                               prefix=badCTFName,
                                                               copyInfo=True)
            for ctfSerie in ctfSeries:
                ctfSerieClon = ctfSerie.clone()
                if CTFSerieStates.UNCHECKED in self.tree.item(ctfSerie.getTsId(),
                                                              'tags'):
                    # Adding the ctfSerie to the good set of ctfTomoSeries
                    outputSetOfgoodCTFTomoSeries.append(ctfSerieClon)
                    outputSetOfgoodCTFTomoSeries.setSetOfTiltSeries(self._inputSetOfTiltSeries)

                else:
                    # Adding the ctfSerie to the bad set of ctfTomoSeries
                    outputSetOfbadCTFTomoSeries.append(ctfSerieClon)
                    outputSetOfbadCTFTomoSeries.setSetOfTiltSeries(self._inputSetOfTiltSeries)

                for item in ctfSerie.iterItems():
                    ctfEstItem = item.clone()
                    ctfSerieClon.append(ctfEstItem)

            outputgoodCTFSetName = 'goodSetOfCTFTomoSeries%s' % suffix
            outputbadCTFSetName = 'badSetOfCTFTomoSeries%s' % suffix

            if len(outputSetOfgoodCTFTomoSeries) > 0:
                protocol._defineOutputs(**{outputgoodCTFSetName: outputSetOfgoodCTFTomoSeries})

            if len(outputSetOfbadCTFTomoSeries) > 0:
                protocol._defineOutputs(**{outputbadCTFSetName: outputSetOfbadCTFTomoSeries})

            protocol._store()
            self.cancel()

    def _showHelp(self, event=None):
        showInfo('CTFTomoSeries viewer help',
                 'This viewer calculates the standard deviation with respect '
                 'to the mean of the defocusU and defocusV values. If the '
                 'values of the images are not in the 20% range from the average '
                 'they are marked as *Failed* and therefore the CTFTomoSerie is '
                 'marked as *Failed* as well.\n\n'
                 'On the other hand, the viewer allows you to create two '
                 'subsets of CTFTomoSeries which are classified as good '
                 'and bad respectively.\n\n'
                 'Note: The series that are checked are the ones that '
                 'represent the bad CTFTomoSeries', self.parent)

    def _getSuffix(self, protocol):
        """
        Return the number of the last output in order to complete the new
        output with a suffix
        """
        maxCounter = -1
        pattern = 'goodSetOfCTFTomoSeries'
        for attrName, _ in protocol.iterOutputAttributes():
            suffix = attrName.replace(pattern, '')
            try:
                counter = int(suffix)
            except:
                counter = 1  # when there is no number, assume 1
            maxCounter = max(counter, maxCounter)

        return str(maxCounter + 1) if maxCounter > 0 else ''

    def _createBottomPanel(self, bottomPanel):
        self._createCTFEstimationGUI(bottomPanel)
        self.initial_focus = self.tree

    def _createCTFEstimationGUI(self, bottomPanel):
        # Create a division Paned
        pw = tk.PanedWindow(bottomPanel, orient=tk.HORIZONTAL)
        # Create a left panel to put the tree
        bottomleftPanel = tk.Frame(pw)
        bottomleftPanel.grid(row=0, column=0, padx=0, pady=0, sticky='news')
        self._createTree(bottomleftPanel)
        pw.add(bottomleftPanel)
        # Panel to put the plotter
        self.bottomRightPanel = ttk.Frame(pw)
        self.bottomRightPanel.grid(row=0, column=1, padx=0, pady=0, sticky='news')
        self._createPloter(self.bottomRightPanel)
        pw.add(self.bottomRightPanel)
        pw.pack(fill=BOTH, expand=True)
        # This method is used to show sash
        pw.configure(sashrelief=RAISED)

    def _createTree(self, parent):

        gui.configureWeigths(parent)

        self.tree = CTFEstimationTree(parent, self.provider,
                                      selectmode=self._selectmode)
        item = self.tree.identify_row(0)
        self.tree.selection_set(item)
        self.tree.focus(item)
        self.tree.selectedItem = item
        self.im_checked = gui.getImage(Icon.CHECKED)
        self.im_unchecked = gui.getImage(Icon.UNCHECKED)
        self.tree.tag_configure(CTFSerieStates.UNCHECKED,
                                image=self.im_unchecked)
        self.tree.tag_configure(CTFSerieStates.CHECKED,
                                image=self.im_checked)
        self.tree.tag_configure(CTFSerieStates.EVEN, background='#F2F2F2',
                                foreground='black')
        self.tree.tag_configure(CTFSerieStates.ODD, background='#E6E6E6',
                                foreground='black')
        self.tree.bind("<Button-1>", self._createPloter, True)

    def _createPloter(self, event):
        itemSelected = self.tree.getSelectedItem()
        obj = self.tree.getSelectedObj()
        self._checkedItems = self.tree._checkedItems
        if self._checkedItems and self._checkedItems != len(self.provider.getCTFSeries()):
            self.generateSubsetButton['state'] = tk.NORMAL
        else:
            self.generateSubsetButton['state'] = tk.DISABLED

        if self.tree.parent(itemSelected):  # child item
            if obj is not None:
                plotterPanel = tk.Frame(self.bottomRightPanel)

                for ctfSerie in self.provider.getCTFSeries():
                    if ctfSerie.getTsId() in itemSelected:
                        ctfId = int(itemSelected.split('.')[-1])
                        # TODO: is ctfSerie ordered by id?
                        fig = self._show2DPLot(ctfSerie, ctfId)
                        canvas = FigureCanvasTkAgg(fig, master=plotterPanel)
                        canvas.draw()
                        canvas.get_tk_widget().pack(fill=BOTH, expand=0)
                        plotterPanel.grid(row=0, column=1, sticky='news')
                        break

        else:  # parent item
            if obj is not None:
                plotterPanel = tk.Frame(self.bottomRightPanel)
                angList = []
                defocusUList = []
                phShList = []
                resList = []

                for ctfSerie in self.provider.getCTFSeries():
                    ts = ctfSerie.getTiltSeries()
                    if ctfSerie.getTsId() == itemSelected:
                        for item in ctfSerie.iterItems(orderBy='id'):
                            index = int(item.getIndex())
                            angList.append(int(ts[index].getTiltAngle()))
                            defocusUList.append(item.getDefocusU())
                            phShList.append(item.getPhaseShift() if item.hasPhaseShift() else 0)
                            resList.append(item.getResolution())

                        fig = Figure(figsize=(7, 7), dpi=100)
                        defocusPlot = fig.add_subplot(111)
                        defocusPlot.grid()
                        defocusPlot.set_title(itemSelected)
                        defocusPlot.set_xlabel('Tilt angle')
                        defocusPlot.set_ylabel('Defocus', color='tab:red')
                        defocusPlot.plot(angList, defocusUList, marker='.',
                                         color='tab:red', label='Defocus (A)')

                        if item.hasPhaseShift():
                            phShPlot = defocusPlot.twinx()
                            phShPlot.set_ylim(0, 180)
                            phShPlot.set_ylabel('Phase shift', color='tab:green')
                            phShPlot.plot(angList, phShList, marker='.',
                                          color='tab:green', label='Phase shift (deg)')
                        else:  # no phase shift, plot resolution
                            resPlot = defocusPlot.twinx()
                            resPlot.set_ylim(0, 30)
                            resPlot.set_ylabel('Resolution', color='tab:green')
                            resPlot.plot(angList, resList, marker='.',
                                         color='tab:green', label='Resolution (A)')

                        fig.legend()
                        canvas = FigureCanvasTkAgg(fig, master=plotterPanel)
                        canvas.draw()
                        canvas.get_tk_widget().pack(fill=BOTH, expand=0)
                        plotterPanel.grid(row=0, column=1, sticky='news')
                        break
