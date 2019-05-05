# PYTHON
import ctypes
import sys
import numpy as np
import scipy.constants as scicon
import time
# OWN ONES
from process_ac import *
# EXTERNAL
from scipy.optimize import minimize
from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5

from PyQt5.QtWinExtras import QWinTaskbarButton
from PyQt5.QtGui import QIcon, QFont, QDoubleValidator
from PyQt5.QtWidgets import (QMainWindow, QWidget, QApplication, QPushButton, QLabel,
                             QDoubleSpinBox, QFormLayout, QCheckBox, QVBoxLayout, QMessageBox,
                             QHBoxLayout, QFileDialog, QDialog, QLineEdit, QListWidget, QListWidgetItem, QTabWidget)
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

class ACGui(QMainWindow):

    def __init__(self):
    
        super().__init__()
        
        self.initUI()
        
    def initUI(self):
        
        """ Things to do with how the window is shown """
        
        self.setWindowTitle('AC Processing')
        self.setWindowIcon(QIcon('double_well_potential_R6p_icon.ico'))
        
        #self.taskbar_btn = QWinTaskbarButton()
        #self.taskbar_btn.setWindow(super)
        
        """ FONT STUFF """
        self.headline_font = QFont()
        self.headline_font.setBold(True)
        
        """ Setting up the main tab widget """
        self.all_the_tabs = QTabWidget()
        self.setCentralWidget(self.all_the_tabs)
        
        """ Constructing the data analysis tab """
        self.data_analysis_tab = QWidget()
        self.all_the_tabs.addTab(self.data_analysis_tab, 'Analysis')
        
        self.main_layout = QHBoxLayout()
        
        # Data containers
        self.data_file_origin = None
        self.data_T = None
        self.data_tau = None
        
        self.data_used_pointer = None
        self.data_not_used_pointer = None
        self.plot_of_fit_pointer = None
        
        self.used_T = None
        self.not_used_T = None
        
        self.used_tau = None
        self.not_used_tau = None
        
        self.fitted_parameters = None
        self.used_indices = None
        
        self.simulation_items = []
        
        """ Adding load controls """
        self.load_layout = QVBoxLayout()
        self.load_layout.addStretch(1)
        
        self.see_fit_btn = QPushButton('Fitted params')
        self.see_fit_btn.clicked.connect(self.print_fitted_params)
        self.load_layout.addWidget(self.see_fit_btn)
        
        self.load_btn = QPushButton('Load')
        self.load_btn.clicked.connect(self.load_t_tau_data)
        self.load_layout.addWidget(self.load_btn)
        
        self.reset_axes_btn = QPushButton('Reset axes')
        self.reset_axes_btn.clicked.connect(self.see_all_on_axes)
        self.load_layout.addWidget(self.reset_axes_btn)
        
        """ Adding plotting controls """
        self.plot_layout = QVBoxLayout()
        
        # Adding the canvas for plotting relaxation times
        self.tau_t_fig = Figure(figsize=(5,3))
        self.tau_t_canvas = FigureCanvas(self.tau_t_fig)
        self.tau_t_tools = NavigationToolbar(self.tau_t_canvas, self)
        self.tau_t_ax = self.tau_t_fig.add_subplot(111)
        
        # Plotting in the axes
        self.tau_t_ax.set_xlabel('Temperature [$K^{-1}$]')
        self.tau_t_ax.set_ylabel(r'$\ln{\tau}$ [$\ln{s}$]')
        
        self.plot_layout.addWidget(self.tau_t_canvas)
        self.plot_layout.addWidget(self.tau_t_tools)
        
        # Adding fit controls
        self.fit_layout = QVBoxLayout()
        
        self.cb_headline = QLabel('Fit types to consider')
        self.cb_headline.setFont(self.headline_font)
        self.fit_layout.addWidget(self.cb_headline)
        
        self.orbach_cb = QCheckBox('Orbach')
        self.orbach_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_layout.addWidget(self.orbach_cb)
        
        self.raman_cb = QCheckBox('Raman')
        self.raman_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_layout.addWidget(self.raman_cb)
        
        self.qt_cb = QCheckBox('QT')
        self.qt_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_layout.addWidget(self.qt_cb)
        
        # Adding temperature controls
        self.temp_headline = QLabel('Temperature interval')
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
        
        # Adding a list to hold information about simulations
        self.simulations_headline = QLabel('Simulations')
        self.simulations_headline.setFont(self.headline_font)
        self.fit_layout.addWidget(self.simulations_headline)
        
        self.list_of_simulations = QListWidget()
        self.list_of_simulations.doubleClicked.connect(self.edit_a_simulation)
        self.fit_layout.addWidget(self.list_of_simulations)
        
        # Adding buttons to control simulation list
        self.sim_btn_layout = QHBoxLayout()
        
        self.delete_sim_btn = QPushButton('Delete')
        self.delete_sim_btn.clicked.connect(self.delete_sim)
        self.sim_btn_layout.addWidget(self.delete_sim_btn)
        
        self.edit_sim_btn = QPushButton('Edit')
        self.edit_sim_btn.clicked.connect(self.edit_a_simulation)
        self.sim_btn_layout.addWidget(self.edit_sim_btn)
        
        self.new_sim_btn = QPushButton('New')
        self.new_sim_btn.clicked.connect(self.add_new_simulation)
        self.sim_btn_layout.addWidget(self.new_sim_btn)
        
        self.fit_layout.addLayout(self.sim_btn_layout)
        
        # Finalizing layout
        self.main_layout.addLayout(self.load_layout)
        self.main_layout.addLayout(self.plot_layout)
        self.main_layout.addLayout(self.fit_layout)
        self.data_analysis_tab.setLayout(self.main_layout)
        
        """ Creating the data treatment tab """
        self.data_treatment_tab = QWidget()
        self.all_the_tabs.addTab(self.data_treatment_tab, 'Data treatment')
        
        self.data_treatment_layout = QHBoxLayout()
        
        self.data_treatment_lbl = QLabel('Here is where the data treatment is supposed to go on!')
        self.data_treatment_layout.addWidget(self.data_treatment_lbl)
        
        self.data_treatment_tab.setLayout(self.data_treatment_layout)
        
        # Showing the GUI
        self.show()
    
    def print_fitted_params(self):
    
        print(self.fitted_parameters)
    
    def see_all_on_axes(self):
        
        s = 0
        if len(self.tau_t_ax.lines)<1: pass
        else:
            while True:
                start = self.tau_t_ax.lines[s]
                if len(start.get_xdata())<1:
                    s += 1
                else:
                    break
        
            x = start.get_xdata()
            y = start.get_ydata()
            
            new_x = [x.min(), x.max()]
            new_y = [y.min(), y.max()]
            
            for i in range(s+1, len(self.tau_t_ax.lines)):
                x = self.tau_t_ax.lines[i].get_xdata()
                y = self.tau_t_ax.lines[i].get_ydata()
                
                if len(x)>1 and len(y)>1:
                    if x.min()<new_x[0]: new_x[0] = x.min()
                    if x.max()>new_x[1]: new_x[1] = x.max()
                    if y.min()<new_y[0]: new_y[0] = y.min()
                    if y.max()>new_y[1]: new_y[1] = y.max()
            
            self.tau_t_ax.set_xlim(new_x[0]-0.1*(new_x[1]-new_x[0]),new_x[1]+0.1*(new_x[1]-new_x[0]))
            self.tau_t_ax.set_ylim(new_y[0]-0.1*(new_y[1]-new_y[0]),new_y[1]+0.1*(new_y[1]-new_y[0]))
            self.tau_t_canvas.draw()
    
    def add_new_simulation(self):
    
        sim_dialog = SimulationDialog(fitted_parameters=self.fitted_parameters,
                                      plot_type_list=[],
                                      plot_parameters={'tQT': 0.01, 'Cr': 0.00, 'n': 0.00, 't0': 0.00, 'Ueff': 0.00},
                                      min_and_max_temps=[0,0])
        
        finished_value = sim_dialog.exec_()
        
        if finished_value:
            
            plot_type = sim_dialog.plot_type_list
            if len(plot_type)<1:
                pass
            else:
                p_fit = sim_dialog.plot_parameters
                T_vals = sim_dialog.min_and_max_temps
                
                plot_to_make = ''.join(plot_type)
                new_item_text = '{}, ({},{}), tQT: {}, Cr: {}, n: {}, t0: {}, Ueff: {}'.format(
                                plot_type, T_vals[0], T_vals[1], p_fit['tQT'], p_fit['Cr'],
                                p_fit['n'], p_fit['t0'], p_fit['Ueff'])
                
                new_list_item = QListWidgetItem()
                
                line = addPartialModel(self.tau_t_fig,
                                        T_vals[0],
                                        T_vals[1],
                                        self.prepare_sim_dict_for_plotting(p_fit),
                                        plotType=plot_to_make)
                                        
                list_item_data = {'plot_type': plot_type,
                                  'p_fit': p_fit,
                                  'T_vals': T_vals,
                                  'line': line}
                
                self.list_of_simulations.addItem(new_list_item)
                new_list_item.setText(new_item_text)
                new_list_item.setData(32, list_item_data)
                
                self.tau_t_canvas.draw()
                
        else:
            pass
    
    def edit_a_simulation(self):
        
        try:
            sim_item = self.list_of_simulations.selectedItems()[0]
        except IndexError:
            pass
        else:
            # Reading off information from the selected item
            old_data = sim_item.data(32)
            old_plot_type_input = old_data['plot_type']
            old_p_fit = old_data['p_fit']
            old_T_vals = old_data['T_vals']
            old_line = old_data['line']
            
            
            
            # Opening simulation dialog with old parameters
            sim_dialog = SimulationDialog(fitted_parameters=self.fitted_parameters,
                                          plot_type_list=old_plot_type_input,
                                          plot_parameters=old_p_fit,
                                          min_and_max_temps=old_T_vals)
            
            finished_value = sim_dialog.exec_()
            
            if finished_value:
                # Reading new parameters of simulation
                new_plot_type = sim_dialog.plot_type_list
                
                if len(new_plot_type)<1:
                    pass
                else:
                    new_p_fit = sim_dialog.plot_parameters
                    new_T_vals = sim_dialog.min_and_max_temps
                    
                    plot_to_make = ''.join(new_plot_type)
                    new_item_text = '{}, ({},{}), tQT: {}, Cr: {}, n: {}, t0: {}, Ueff: {}'.format(
                                    new_plot_type, new_T_vals[0], new_T_vals[1], new_p_fit['tQT'], new_p_fit['Cr'],
                                    new_p_fit['n'], new_p_fit['t0'], new_p_fit['Ueff'])
                    
                    self.tau_t_ax.lines.remove(old_line)
                    
                    new_line = addPartialModel(self.tau_t_fig,
                                               new_T_vals[0],
                                               new_T_vals[1],
                                               self.prepare_sim_dict_for_plotting(new_p_fit),
                                               plotType=plot_to_make)
                    
                    list_item_data = {'plot_type': new_plot_type,
                                      'p_fit': new_p_fit,
                                      'T_vals': new_T_vals,
                                      'line': new_line}
                    
                    self.tau_t_canvas.draw()
                    
                    sim_item.setData(32, list_item_data)
                    sim_item.setText(new_item_text)
            else:
                pass
        
    def delete_sim(self):
        
        try:
            sim_item = self.list_of_simulations.selectedItems()[0]
        except IndexError:
            pass
        else:
            line_pointer = sim_item.data(32)['line']
            self.tau_t_ax.lines.remove(line_pointer)
            self.tau_t_canvas.draw()
            
            item_row = self.list_of_simulations.row(sim_item)
            sim_item = self.list_of_simulations.takeItem(item_row)
            
            del sim_item
        
    def plot_t_tau_on_axes(self):
        
        if self.data_used_pointer is not None:
            self.tau_t_ax.lines.remove(self.data_used_pointer)
        if self.data_not_used_pointer is not None:
            self.tau_t_ax.lines.remove(self.data_not_used_pointer)
        
        self.data_used_pointer, = self.tau_t_ax.plot(1/self.used_T, np.log(self.used_tau), 'bo')
        self.data_not_used_pointer, = self.tau_t_ax.plot(1/self.not_used_T, np.log(self.not_used_tau), 'ro')
        
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
            print('Encountered value error... check your input file!')
        except OSError:
            print('File was not found! Empty file name?')
        else:
            self.read_indices_for_used_temps()
            self.plot_t_tau_on_axes()
    
    def prepare_sim_dict_for_plotting(self, p_fit_gui_struct):

        params = []
        quantities = []
        sigmas = [0]*5

        for key, val in p_fit_gui_struct.items():
            params.append(val)
            quantities.append(key)
        
        Ueff = params[quantities.index('Ueff')]
        params[quantities.index('Ueff')] = Ueff*scicon.Boltzmann
        
        p_fit_script_type = {'params': params,
                             'quantities': quantities,
                             'sigmas': sigmas}
        
        return p_fit_script_type
    
    def read_fit_type_cbs(self):
    
        list_of_checked = []
        if self.qt_cb.isChecked(): list_of_checked.append('QT')
        if self.raman_cb.isChecked(): list_of_checked.append('R')
        if self.orbach_cb.isChecked(): list_of_checked.append('O')
        fitToMake = ''.join(list_of_checked)
        
        return fitToMake
    
    def read_indices_for_used_temps(self):
        
        min_t = self.temp_line[1].value()
        max_t = self.temp_line[3].value()
        
        try:
            self.used_indices = [list(self.data_T).index(t) for t in self.data_T if t>=min_t and t<=max_t]
            
            self.used_T = self.data_T[self.used_indices]
            self.used_tau = self.data_tau[self.used_indices]
            
            self.not_used_T = np.delete(self.data_T, self.used_indices)
            self.not_used_tau = np.delete(self.data_tau, self.used_indices)
            
        except (AttributeError, TypeError):
            print('No data have been selected yet!')
        
    
    def set_new_temp_ranges(self):
    
        new_max_for_low = self.temp_line[3].value()
        new_min_for_high = self.temp_line[1].value()
        self.temp_line[1].setRange(0,new_max_for_low)
        self.temp_line[3].setRange(new_min_for_high,1000)
        
        self.read_indices_for_used_temps()
        if self.data_T is not None:
            self.plot_t_tau_on_axes()
        
    def make_the_fit(self):
        
        try:
            Tmin = self.temp_line[1].value()
            Tmax = self.temp_line[3].value()
            perform_this_fit = self.read_fit_type_cbs()
            assert Tmin != Tmax
            assert perform_this_fit != ''
            
            fig3, p_fit = fitRelaxation(self.data_file_origin, (Tmin, Tmax), fitType=perform_this_fit)
            self.fitted_parameters = p_fit
        
        except (AssertionError, RuntimeError, ValueError) as error:
            
            error_type = error.__class__.__name__
            if error_type == 'AssertionError':
                msg_text = 'Check your temperature and fit settings'
                msg_details = """Possible errors:
 - min and max temperatures are the same
 - no fit options have been selected"""
            elif error_type == 'RuntimeError':
                msg_text = 'This fit cannot be made within the set temperatures'
                msg_details = ''
            elif error_type == 'ValueError':
                msg_text = 'No file has been loaded'
                msg_details = ''
            
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle('Fit aborted')
            msg.setText(msg_text)
            msg.setDetailedText(msg_details)
            msg.exec_()

class SimulationDialog(QDialog):

    def __init__(self,
                 parent=None,
                 fitted_parameters=None,
                 plot_type_list=[],
                 plot_parameters={'tQT': 0.01, 'Cr': 0.00, 'n': 0.00, 't0': 0.00, 'Ueff': 0.00},
                 min_and_max_temps=[0]*2):
    
        super(SimulationDialog, self).__init__()
        
        self.setWindowTitle('Add simulation')
        
        # Containers for objects
        self.lineedit_inputs = {}
        
        self.headline_font = QFont()
        self.headline_font.setBold(True)
        
        self.plot_type_list = plot_type_list
        self.plot_parameters = plot_parameters
        
        self.min_and_max_temps = min_and_max_temps
        
        self.layout = QVBoxLayout()
        
        if fitted_parameters is not None:
            self.use_fit_cb = QCheckBox('Use fitted parameters')
            self.use_fit_cb.stateChanged.connect(self.fit_take_control)
            self.layout.addWidget(self.use_fit_cb)
            
            self.fitted_parameters = {}
            for quant in fitted_parameters['quantities']:
                try:
                    idx = fitted_parameters['quantities'].index(quant)
                    param = fitted_parameters['params'][idx]
                except IndexError:
                    param = None
                finally:
                    if quant=='Ueff':
                        param/=scicon.Boltzmann
                    self.fitted_parameters[quant] = param
            
        # Controls to play with temperature
        self.temp_headline = QLabel('Temperature')
        self.temp_headline.setFont(self.headline_font)
        self.layout.addWidget(self.temp_headline)
        
        self.temp_hbl = QHBoxLayout()
        
        self.temp_min = QDoubleSpinBox()
        self.temp_min.setValue(min_and_max_temps[0])
        self.temp_min.editingFinished.connect(self.temp_interval_changed)
        self.temp_hbl.addWidget(self.temp_min)
        
        self.temp_max = QDoubleSpinBox()
        self.temp_max.setValue(min_and_max_temps[1])
        self.temp_max.editingFinished.connect(self.temp_interval_changed)
        self.temp_hbl.addWidget(self.temp_max)
        
        self.layout.addLayout(self.temp_hbl)
        
        # Controls for which type of plot to consider
        self.plot_headline = QLabel('Plot type to make')
        self.plot_headline.setFont(self.headline_font)
        self.layout.addWidget(self.plot_headline)
        
        self.plot_type_hbl = QHBoxLayout()
        
        self.use_orbach = QCheckBox('Orbach')
        if 'O' in self.plot_type_list: self.use_orbach.setChecked(True)
        self.use_orbach.clicked.connect(self.plot_type_changed)
        self.plot_type_hbl.addWidget(self.use_orbach)
        
        self.use_raman = QCheckBox('Raman')
        if 'R' in self.plot_type_list: self.use_raman.setChecked(True)
        self.use_raman.clicked.connect(self.plot_type_changed)
        self.plot_type_hbl.addWidget(self.use_raman)
        
        self.use_qt = QCheckBox('QT')
        if 'QT' in self.plot_type_list: self.use_qt.setChecked(True)
        self.use_qt.clicked.connect(self.plot_type_changed)
        self.plot_type_hbl.addWidget(self.use_qt)
        
        self.layout.addLayout(self.plot_type_hbl)
        
        self.sim_vals_layout = QFormLayout()
        
        self.tqt_val = QLineEdit()
        self.lineedit_inputs['tQT'] = self.tqt_val
        self.tqt_val.setValidator(QDoubleValidator())
        self.tqt_val.setText(str(self.plot_parameters['tQT']))
        self.tqt_val.editingFinished.connect(self.param_values_changed)
        self.sim_vals_layout.addRow('t_QT', self.tqt_val)
        
        self.cr_val = QLineEdit()
        self.lineedit_inputs['Cr'] = self.cr_val
        self.cr_val.setValidator(QDoubleValidator())
        self.cr_val.setText(str(self.plot_parameters['Cr']))
        self.cr_val.editingFinished.connect(self.param_values_changed)
        self.sim_vals_layout.addRow('C_R', self.cr_val)
        
        self.n_val = QLineEdit()
        self.lineedit_inputs['n'] = self.n_val
        self.n_val.setValidator(QDoubleValidator())
        self.n_val.setText(str(self.plot_parameters['n']))
        self.n_val.editingFinished.connect(self.param_values_changed)
        self.sim_vals_layout.addRow('n', self.n_val)
        
        self.t0_val = QLineEdit()
        self.lineedit_inputs['t0'] = self.t0_val
        self.t0_val.setValidator(QDoubleValidator())
        self.t0_val.setText(str(self.plot_parameters['t0']))
        self.t0_val.editingFinished.connect(self.param_values_changed)
        self.sim_vals_layout.addRow('t0', self.t0_val)
        
        self.Ueff_val = QLineEdit()
        self.lineedit_inputs['Ueff'] = self.Ueff_val
        self.Ueff_val.setValidator(QDoubleValidator())
        self.Ueff_val.setText(str(self.plot_parameters['Ueff']))
        self.Ueff_val.editingFinished.connect(self.param_values_changed)
        self.sim_vals_layout.addRow('U_eff', self.Ueff_val)
        
        self.layout.addLayout(self.sim_vals_layout)
        
        # Making control buttons at the end
        self.button_layout = QHBoxLayout()
        
        self.cancel_btn = QPushButton('Cancel')
        self.cancel_btn.setAutoDefault(False)
        self.cancel_btn.clicked.connect(self.reject)
        self.button_layout.addWidget(self.cancel_btn)
        
        self.accept_btn = QPushButton('Ok')
        self.accept_btn.clicked.connect(self.accept)
        self.button_layout.addWidget(self.accept_btn)
        
        self.layout.addLayout(self.button_layout)
        
        self.setLayout(self.layout)
        
        self.show()
    
    def fit_take_control(self):
        
        if self.use_fit_cb.isChecked():
            
            for key, val in self.fitted_parameters.items():
                if val is not None:
                    self.lineedit_inputs[key].setText('{}'.format(val))
                self.lineedit_inputs[key].setReadOnly(True)
            
        else:
            self.tqt_val.setReadOnly(False)
            self.cr_val.setReadOnly(False)
            self.n_val.setReadOnly(False)
            self.t0_val.setReadOnly(False)
            self.Ueff_val.setReadOnly(False)
        
        self.param_values_changed()
    
    def param_values_changed(self):
        
        self.plot_parameters['tQT'] = float(self.tqt_val.text())
        self.plot_parameters['Cr'] = float(self.cr_val.text())
        self.plot_parameters['n'] = float(self.n_val.text())
        self.plot_parameters['t0'] = float(self.t0_val.text())
        self.plot_parameters['Ueff'] = float(self.Ueff_val.text())
        
    def plot_type_changed(self):
        
        self.plot_type_list = []
        if self.use_qt.isChecked(): self.plot_type_list.append('QT')
        if self.use_raman.isChecked(): self.plot_type_list.append('R')
        if self.use_orbach.isChecked(): self.plot_type_list.append('O')
        
    def temp_interval_changed(self):
        
        try:
            self.min_and_max_temps[0] = self.temp_min.value()
            self.min_and_max_temps[1] = self.temp_max.value()
            assert self.min_and_max_temps[0]<=self.min_and_max_temps[1]
        except AssertionError:
            pass
            
if __name__ == '__main__':
    
    myappid = 'AC Processing v1.0'
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    
    app = QApplication(sys.argv)
    w = ACGui()
    sys.exit(app.exec_())
        
        
        