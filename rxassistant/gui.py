import sys
import matplotlib
# matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import numpy as np

from refnx.reflect import SLD, ReflectModel


class RxCanvas(FigureCanvasQTAgg):
    wavelength = 1.542
    coll_dist = 1200

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(RxCanvas, self).__init__(fig)
        self.axes.set_xlabel("$Q$ [$\\AA^{-1}$]")
        self.axes.set_ylabel("$R$")
        self.axes.set_yscale("log")
        self.line = None

    def plot(self, theta, dtheta=None, flux=None, acq_time=None):
        if flux is None or flux == 0:
            flux = 0
        if dtheta is None:
            dtheta = 0
        if acq_time is not None:
            acq_time = acq_time[theta > 0]
        else:
            acq_time = 0
        theta = theta[theta > 0]
        if len(theta) == 0:
            return
        si = SLD(20.062 - 1j * 0.458, name="Si")
        air = SLD(0, name="air")
        si_layer = si(0, 0)
        air_layer = air(0, 0)
        structure = air_layer | si_layer
        q = 4 * np.pi * np.sin(np.deg2rad(theta)) / self.wavelength
        dq = q *dtheta / theta
        model = ReflectModel(structure, dq=0, bkg=0)
        r = model(q, x_err=dq)
        dr = r / np.sqrt(r*flux*acq_time)
        dr[np.isinf(dr)] = 0

        if self.line is None:
            self.line = self.axes.errorbar(q, r, yerr=dr, fmt='.')
        else:
            self.line.remove()
            self.line = self.axes.errorbar(q, r, yerr=dr, fmt='b.')
            # self.line.set_data(theta, r)
        self.axes.set_xlim(0, np.max(q))
        self.axes.set_ylim(np.min(r), 2)
        self.draw()


class ParametersWidget(QtWidgets.QGroupBox):
    parametersChanged = QtCore.pyqtSignal(dict)

    def __init__(self, group_name, names=[], types=[]):
        """
        Create ParametersWidget composed of  form layou with labels and QLineEdit
        :param names: parameters names
        :rtype list
        :param types: parameters type (float, int or str)
        :rtype list
        """
        super().__init__(group_name)
        if len(types) != len(names):
            raise IndexError("parameters names list has not the same length as types list")
        self._names = names
        self._types = types
        self.lineEditList = []
        # float validator with "." as decimal delimiter
        float_validator = QtGui.QDoubleValidator()
        locale = QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates)
        float_validator.setLocale(locale)
        layout = QtWidgets.QFormLayout()
        for pname, ptype in zip(names, types):
            widget = QtWidgets.QLineEdit(parent=self)
            if ptype is float:
                widget.setValidator(float_validator)
            elif ptype is int:
                widget.setValidator(QtGui.QIntValidator())
            else:
                pass
            layout.addRow(pname, widget)
        self.setLayout(layout)
        # connect signals
        for i in range(layout.rowCount()):
            widget = layout.itemAt(i, 1).widget()
            widget.textChanged.connect(self.onParametersChanged)

    def onParametersChanged(self):
        ans = self.get_parameters()
        self.parametersChanged.emit(ans)

    def get_parameters(self):
        """
        :return: the typed parameter values as a dictionnary. The keys are the label of the form layout
        :rtype: dict
        """
        pdict = {}
        layout = self.layout()
        for i in range(layout.rowCount()):

            label_widget = layout.itemAt(i, 0).widget()
            edit_widget = layout.itemAt(i, 1).widget()
            p = edit_widget.text()
            if self._types[i] is float:
                try:
                    p = float(p)
                except ValueError:
                    p = 0.0
            elif self._types[i] is int:
                try:
                    p = int(p)
                except ValueError:
                    p = 0
            else:
                pass
            pdict[label_widget.text()] = p
        return pdict


class SampleParametersWidget(ParametersWidget):

    def __init__(self):
        super().__init__("sample", names=["name", "x", "z", "offset"],
                         types=[str, float, float, float])


class AcquisitionParametersWidget(ParametersWidget):

    def __init__(self):
        super().__init__("acquisition", names=["theta0", "N1", "dtheta1", "t1", "N2", "dtheta2", "t2"],
                         types=[float, int, float, float, int, float, float])

    def get_theta(self):
        """
        :return: the array of incident angles corresponding to the settings
        :rtype: 1D array
        """
        p = self.get_parameters()
        try:
            theta1 = np.linspace(p["theta0"], p["theta0"]+p["N1"]*p["dtheta1"], p["N1"])
            theta2 = np.linspace(theta1[-1] + p["dtheta2"], theta1[-1] + p["dtheta2"] + p["N2"] * p["dtheta2"], p["N2"])
            theta = np.concatenate((theta1, theta2))
        except IndexError:
            theta = np.array([0])
        return theta

    def get_acq_time(self):
        """
        :return: the array of acquisition times corresponding to the settings
        :rtype: 1D array
        """
        p = self.get_parameters()
        t1 = p['t1']*np.ones(p['N1'])
        t2 = p['t2'] * np.ones(p['N2'])
        t = np.concatenate((t1, t2))
        if len(t) ==0:
            return None
        else:
            return t

class ConfigurationParameters(ParametersWidget):

    def __init__(self):
        super().__init__("config", names=["flux", "dtheta"], types=[float, float])



class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Create the maptlotlib FigureCanvas object,
        # which defines a single set of axes as self.axes.
        # self.parameters = SampleParametersWidget()
        central_widget = QtWidgets.QWidget()
        self.acq_parameters = AcquisitionParametersWidget()
        self.sample_parameters = SampleParametersWidget()
        self.config_parameter = ConfigurationParameters()
        # self.parameters = ParametersWidget(names=["trucs"], types=[str])

        self.canvas = RxCanvas(self, width=5, height=4, dpi=100)
        # self.setCentralWidget(sc)
        layout = QtWidgets.QHBoxLayout()
        vlayout = QtWidgets.QVBoxLayout()
        vlayout.addWidget(self.sample_parameters)
        vlayout.addWidget(self.acq_parameters)
        vlayout.addWidget(self.config_parameter)

        layout.addLayout(vlayout)
        layout.addWidget(self.canvas)
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)
        # connect signal
        self.acq_parameters.parametersChanged.connect(self.update_plot)
        self.config_parameter.parametersChanged.connect(self.update_plot)
        # init plot
        self.acq_parameters.onParametersChanged()

    def update_plot(self):
        theta = self.acq_parameters.get_theta()
        pconf = self.config_parameter.get_parameters()
        acq_time = self.acq_parameters.get_acq_time()
        print(acq_time)
        self.canvas.plot(theta, flux=pconf['flux'], dtheta=pconf['dtheta'], acq_time=acq_time)


class AcquisiotnTable(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()
        self.table = QtWidgets.QTableWidget()
        self.table.setColumnCount(3)
        # self.table.setRowCount(2)
        self.table.setHorizontalHeaderLabels(["N", "step (°)", "t (s)"])
        self.addStep()
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.table)
        self.setLayout(layout)

    def addStep(self):
        n = self.table.rowCount()
        print(n)
        self.table.insertRow(n)
        edit = QtWidgets.QLineEdit("100")
        edit.setValidator(QtGui.QIntValidator())
        self.table.setCellWidget(n,0,edit)
        edit = QtWidgets.QLineEdit("0.01")
        edit.setValidator(QtGui.QDoubleValidator())
        self.table.setCellWidget(n, 1, edit)
        edit = QtWidgets.QLineEdit("1")
        edit.setValidator(QtGui.QDoubleValidator())
        self.table.setCellWidget(n, 2, edit)
        self.table.setItem(n,0, QtWidgets.QTableWidgetItem("100"))
        self.table.setItem(n,1, QtWidgets.QTableWidgetItem("0.01"))
        self.table.setItem(n,2, QtWidgets.QTableWidgetItem("1"))


    def get_parameters(self):
        nList = []
        dthetaList = []
        tList = []
        for i in range(self.table.rowCount()):
            text = self.table.item(i, 0).text()
            nList.append(int(self.table.item(i, 0).text()))
            dthetaList.append(float(self.table.item(i, 1).text()))
            tList.append()


class AcquisitionWidget(QtWidgets.QWidget):
    parametersChanged = QtCore.pyqtSignal(list)

    def __init__(self):
        super().__init__()
        self.table = QtWidgets.QTableView()
        self.model = AcquisitionTableModel()
        self.table.setModel(self.model)
        # add and remove step button
        self.addStepButton = QtWidgets.QPushButton("+")
        self.removeStepButton = QtWidgets.QPushButton("-")
        # total estimated acquisition time
        self.totalTimeLineEdit = QtWidgets.QLineEdit("")
        self.totalTimeLineEdit.setReadOnly(True)
        # first incident angle
        self.thetaStartLineEdit = QtWidgets.QLineEdit("0")
        self.thetaStartLineEdit.setValidator(QtGui.QDoubleValidator())
        # first incident angle
        self.thetaEndLineEdit = QtWidgets.QLineEdit()
        self.thetaEndLineEdit.setReadOnly(True)
        # layout
        mainLayout = QtWidgets.QVBoxLayout()

        btnlayout = QtWidgets.QHBoxLayout()
        btnlayout.addWidget(self.addStepButton)
        btnlayout.addStretch(1)
        btnlayout.addWidget(self.removeStepButton)
        btnlayout.addStretch(1)

        mainLayout.addLayout(btnlayout)
        mainLayout.addWidget(self.table)
        mainLayout.addWidget(self.totalTimeLineEdit)
        mainLayout.addWidget(self.thetaEndLineEdit)
        mainLayout.addStretch(1)

        self.setLayout(mainLayout)

        # connect signals
        self.addStepButton.clicked.connect(self.onAddClicked)
        self.removeStepButton.clicked.connect(self.onRemoveClicked)
        self.model.dataChanged.connect(self.onParametersChanged)

    def onAddClicked(self):
        n = self.model.rowCount()
        self.model.insertRow(n)

    def onRemoveClicked(self):
        n = self.model.rowCount()
        if n > 0:
            self.model.removeRow(n-1)

    def onParametersChanged(self):
        print("changing")
        t = self.model.get_acq_time()
        try:
            thetaStart = float(self.thetaStartLineEdit.text())
        except ValueError:
            thetaStart = 0.
        theta = self.model.get_theta() + thetaStart
        self.parametersChanged.emit([theta, t])
        self.totalTimeLineEdit.setText(str(np.sum(t)))
        self.thetaEndLineEdit.setText(str(np.max(theta)))



class AcquisitionTableModel(QtCore.QAbstractTableModel):
    def __init__(self):
        super().__init__()
        self._data = [[100, 0.01, 1]]
        self._headers = ["intervals", "step (°)", "t (s)"]

    def rowCount(self, parent=QtCore.QModelIndex()):
        # The length of the outer list.
        return len(self._data)

    def columnCount(self, parent=QtCore.QModelIndex()):
        # The following takes the first sub-list, and returns
        # the length (only works if all rows are an equal length)
        return 3

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.DisplayRole):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            # return f"Column {section + 1}"
            return self._headers[section]
        if orientation == Qt.Vertical and role == Qt.DisplayRole:
            return f"{section+1}"

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole or role == Qt.EditRole:
                value = self._data[index.row()][index.column()]
                return str(value)

    def setData(self, index, value, role=Qt.EditRole):
        if role == Qt.EditRole:
            col = index.column()
            if col == 0:
                try:
                    print(value)
                    self._data[index.row()][index.column()] = int(float(value))
                    self.dataChanged.emit(index, index)
                    return True
                except ValueError:
                    return False
            if col == 1 or col ==2:
                try:
                    self._data[index.row()][index.column()] = float(value)
                    self.dataChanged.emit(index, index)
                    return True
                except ValueError:
                    return False
        return False

    def flags(self, index):
        return Qt.ItemIsSelectable | Qt.ItemIsEnabled | Qt.ItemIsEditable

    def insertRow(self, row, parent=QtCore.QModelIndex()):
        self.beginInsertRows(QtCore.QModelIndex(), row, row)
        self._data.insert(row, [100, 0.01, 1])
        self.endInsertRows()
        return True

    def removeRow(self, row: int, parent: QtCore.QModelIndex = QtCore.QModelIndex()) -> bool:
        self.beginRemoveRows(QtCore.QModelIndex(), row, row)
        del self._data[row]
        self.endRemoveRows()
        return True

    def get_theta(self):
        theta = []
        for i in range(len(self._data)):
            theta.append(np.linspace(0, self._data[i][0]*self._data[i][1], self._data[i][0]))

        theta = np.concatenate(theta)
        return theta

    def get_acq_time(self):
        t = []
        for i in range(len(self._data)):
            t.append(np.ones(self._data[i][0])*self._data[i][2])
        return np.concatenate(t)



if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    # w = MainWindow()
    tab = AcquisitionWidget()
    # tab.setModel(mod)
    # mod.insertRow(0)
    # print(mod.get_acq_time())
    # print(mod.get_theta())

    # w = AcquisiotnTable()
    # w.show()
    tab.show()
    app.exec_()