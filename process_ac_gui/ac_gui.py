# PYTHON
import sys
import numpy as np
import time
# OWN ONES
from process_ac import *
# EXTERNAL
from scipy.optimize import minimize
from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5

from PyQt5.QtGui import QIcon, QFont
from PyQt5.QtWidgets import (QWidget, QApplication, QPushButton, QGridLayout, QLabel, QComboBox, QStackedWidget,
                             QDoubleSpinBox, QFormLayout, QCheckBox, QSpinBox, QVBoxLayout, QMessageBox,
                             QHBoxLayout, QFileDialog)
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    
class ACGui(QWidget):

    def __init__(self):
    
        super().__init__()
        
        self.initUI()
        
    def initUI(self):
        
        self.setWindowTitle('AC processing')
        self.headline_font = QFont()
        self.headline_font.setBold(True)
        
        # Data containers
        self.data_file_origin = None
        self.data_T = None
        self.data_tau = None
        self.fitted_parameters = None
        
        self.grid_layout = QGridLayout()
        
        # Adding load controls
        self.load_layout = QVBoxLayout()
        
        self.load_btn = QPushButton('Load')
        self.load_btn.clicked.connect(self.load_t_tau_data)
        self.load_layout.addWidget(self.load_btn)
        
        # Adding the canvas for plotting relaxation times
        self.tau_t_fig = Figure(figsize=(5,3))
        self.tau_t_canvas = FigureCanvas(self.tau_t_fig)
        self.tau_t_tools = NavigationToolbar(self.tau_t_canvas, self)
        self.tau_t_ax = self.tau_t_fig.add_subplot(111)
        
        # Plotting in the axes
        self.tau_t_ax.set_xlabel('Temperature')
        self.tau_t_ax.set_ylabel(r'$\tau$')
        
        self.grid_layout.addWidget(self.tau_t_canvas,0,1,8,8)
        self.grid_layout.addWidget(self.tau_t_tools,8,1,1,8)
        
        # Adding fit controls
        self.fit_layout = QVBoxLayout()
        
        self.cb_headline = QLabel('Fit types to consider')
        self.cb_headline.setFont(self.headline_font)
        self.fit_layout.addWidget(self.cb_headline)
        
        self.orbach_cb = QCheckBox('Orbach')
        self.orbach_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_layout.addWidget(self.orbach_cb)
        
        self.qt_cb = QCheckBox('QT')
        self.qt_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_layout.addWidget(self.qt_cb)
        
        self.raman_cb = QCheckBox('Raman')
        self.raman_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_layout.addWidget(self.raman_cb)
        
        # Adding temperature controls
        self.temp_headline = QLabel('Temperature')
        self.temp_headline.setFont(self.headline_font)
        self.fit_layout.addWidget(self.temp_headline)
        
        self.temp_horizontal_layout = QHBoxLayout()
        self.temp_line = [QLabel('('), QDoubleSpinBox(), QLabel(','), QDoubleSpinBox(), QLabel(')')]
        
        self.temp_line[1].setRange(0,self.temp_line[3].value())
        self.temp_line[1].setSingleStep(0.1)
        self.temp_line[3].setRange(self.temp_line[1].value(),1000)
        self.temp_line[3].setSingleStep(0.1)
        
        self.temp_line[1].editingFinished.connect(self.set_new_temp_ranges)
        self.temp_line[3].editingFinished.connect(self.set_new_temp_ranges)
        for w in self.temp_line:
            self.temp_horizontal_layout.addWidget(w)
        
        self.fit_layout.addLayout(self.temp_horizontal_layout)
        
        # Adding a button to run a fit
        self.run_fit_btn = QPushButton('Run fit!')
        self.run_fit_btn.clicked.connect(self.make_the_fit)
        self.fit_layout.addWidget(self.run_fit_btn)
        
        # Finalizing layout and showing the GUI
        self.grid_layout.addLayout(self.load_layout,1,0,1,1)
        self.grid_layout.addLayout(self.fit_layout,1,9,1,1)
        self.setLayout(self.grid_layout)
        self.show()
    
    def plot_t_tau_on_axes(self):
        
        self.tau_t_ax.clear()
        self.tau_t_ax.plot(1/self.data_T, np.log(self.data_tau), 'bo')
        self.tau_t_canvas.draw()

    
    def load_t_tau_data(self):
        
        starting_directory = os.getcwd()
        filename = QFileDialog().getOpenFileName(self, 'Open file', starting_directory)
        
        self.data_file_origin = filename[0]
        
        try:
            D = np.loadtxt(self.data_file_origin, skiprows=1)
            self.data_T = D[:,0]
            self.data_tau = D[:,1]
        
        except ValueError:
            print('Encountered valueerror... check your input file!')
        except OSError:
            print('File was not found! Empty file name?')
        else:
            self.plot_t_tau_on_axes()
    
    def read_fit_type_cbs(self):
    
        list_of_checked = []
        if self.qt_cb.isChecked(): list_of_checked.append('QT')
        if self.raman_cb.isChecked(): list_of_checked.append('R')
        if self.orbach_cb.isChecked(): list_of_checked.append('O')
        fitToMake = ''.join(list_of_checked)
        
        return fitToMake
    
    def set_new_temp_ranges(self):
    
        new_max_for_low = self.temp_line[3].value()
        new_min_for_high = self.temp_line[1].value()
        self.temp_line[1].setRange(0,new_max_for_low)
        self.temp_line[3].setRange(new_min_for_high,1000)
        
        fit = self.read_fit_type_cbs()
        
    def make_the_fit(self):
        
        try:
            Tmin = self.temp_line[1].value()
            Tmax = self.temp_line[3].value()
            perform_this_fit = self.read_fit_type_cbs()
            assert Tmin != Tmax
            assert perform_this_fit != ''
        
        except AssertionError:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle('Fit aborted')
            msg.setText('Check your temperature and fit settings')
            msg.setDetailedText("""Possible errors:
 - min and max temperatures are the same
 - no fit options have been selected""")
            msg.exec_()
        
        else:
            fig3, p_fit = fitRelaxation('Dy_DBM_ac0Oe.dat', (Tmin, Tmax), fitType=perform_this_fit)
            self.fitted_parameters = p_fit
        
    
    def printSomething(self):
    
        print('Connection established!')
        
if __name__ == '__main__':

    app = QApplication(sys.argv)
    w = ACGui()
    sys.exit(app.exec_())
        
        
        