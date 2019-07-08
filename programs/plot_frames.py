import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os

from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5

from PyQt5.QtCore import Qt, QThread, Signal, QObject
from PyQt5.QtWidgets import (QWidget, QApplication, QPushButton, QGridLayout, QLabel, QComboBox, QStackedWidget, QSlider,
                             QDoubleSpinBox, QFormLayout, QCheckBox, QSpinBox, QVBoxLayout, QMessageBox, QFileDialog,
                             QMainWindow, QHBoxLayout, QVBoxLayout, QTabWidget, QProgressBar)
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

class SixT2Visual(QMainWindow):
    
    ready_to_plot_again = Signal()
    
    def __init__(self):
    
        super().__init__()
        
        self.initUI()
        
    def initUI(self):
        
        # Variables
        self.N = 1
        self.nframes = 0
        self.filepath = None
        self.last_opened_dir = os.getcwd()
        self.frame_x = 0
        self.frame_y = 0
        
        # Main window
        self.setWindowTitle('Data visualization')
        
        self.tab_control = QTabWidget()
        self.setCentralWidget(self.tab_control)
        
        """Design of raw frame tab"""
        self.raw = QWidget()
        self.raw_lo = QVBoxLayout()
        self.tab_control.addTab(self.raw, 'Frames')
        
        # Button controls
        self.btn_lo = QHBoxLayout()
        self.btn_lo.addStretch()
        
        self.left_btn = QPushButton('<')
        self.left_btn.clicked.connect(self.N_subtract)
        self.btn_lo.addWidget(self.left_btn)
        
        self.frame_lbl = QLabel()
        self.btn_lo.addWidget(self.frame_lbl)
        
        self.right_btn = QPushButton('>')
        self.right_btn.clicked.connect(self.N_add)
        self.btn_lo.addWidget(self.right_btn)
        
        self.btn_lo.addStretch()
        
        # Video controls
        self.video_lo = QHBoxLayout()
        
        
        self.video_lo.addStretch()
        self.start_btn = QPushButton('Start')
        self.start_btn.clicked.connect(self.start_frame_video)
        self.video_lo.addWidget(self.start_btn)
        
        self.end_btn = QPushButton('End')
        self.end_btn.clicked.connect(self.stop_frame_video)
        self.video_lo.addWidget(self.end_btn)
        self.video_lo.addStretch()
        
        # Plotting controls
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.axD = self.fig.add_subplot(121)
        self.axU = self.fig.add_subplot(122)
        
        self.rsD = RectangleSelector(self.axD, self.show_sum_in_rectangle,
                                     drawtype='box', useblit=True,
                                     button=1,
                                     minspanx=5, minspany=5,
                                     spancoords='pixels',
                                     interactive=True)
        self.canvas.mpl_connect('key_press_event', self.rsD)
        
        # Meta-functions
        self.data_lo = QHBoxLayout()
        
        self.load_xml_btn = QPushButton('Load')
        self.load_xml_btn.clicked.connect(self.find_file)
        self.data_lo.addWidget(self.load_xml_btn)
        
        self.data_lo.addStretch()
        
        self.progress_bar = QProgressBar()
        self.data_lo.addWidget(self.progress_bar)
        
        # Setting layout
        self.raw_lo.addLayout(self.btn_lo)
        self.raw_lo.addLayout(self.video_lo)
        self.raw_lo.addWidget(self.canvas)
        self.raw_lo.addWidget(self.toolbar)
        self.raw_lo.addLayout(self.data_lo)
        
        # Initialization of raw
        self.raw.setLayout(self.raw_lo)
        self.show()
    
    def show_sum_in_rectangle(self, eclick, erelease):
    
        window = 15
    
        xmin, xmax = int(eclick.xdata), int(np.ceil(erelease.xdata))
        ymin, ymax = int(eclick.ydata), int(np.ceil(erelease.ydata))
        
        #min_sub, max_add = window, window
        #if self.N-window<0: min_sub = self.N
        #if self.N+window>self.nframes: max_add = self.nframes-self.N
        #
        #temp = [{d.attrib['name']:
        #         np.fromstring(d.text, sep=';').reshape(128,128)
        #         for d in f
        #         if d.tag=='Data'}
        #         for f in self.frames[self.N-min_sub:self.N+max_add]]        
        #
        #yD = [np.sum(np.array(t['IntensityDown'])[ymin:ymax,xmin:xmax]) for t in temp]
        #yU = [np.sum(np.array(t['IntensityUp'])[ymin:ymax,xmin:xmax]) for t in temp]
        
        #x = list(range(self.N-min_sub,self.N+max_add))
        #plt.plot(x, yD, label='Down')
        #plt.plot(x, yU, label='Up')
        
        plt.plot([np.sum(d['IntensityUp'][ymin:ymax, xmin:xmax]) for d in self.frame_data])
        plt.plot([np.sum(d['IntensityDown'][ymin:ymax, xmin:xmax]) for d in self.frame_data])
        plt.show()
    
    def find_file(self):
    
        w = QFileDialog()
        w.setDirectory(self.last_opened_dir)
        self.filepath = w.getOpenFileName()[0]
        
        path, ext = os.path.splitext(self.filepath)
        self.last_opened_dir = path
        
        if ext == '.xml':
            self.readXML()
        else:
            pass
        
    def update_N_lbl(self):
    
        self.frame_lbl.setText('{}/{}'.format(self.N, self.nframes))
    
    def N_add(self):
        
        if self.N == self.nframes:
            pass
        else:
            self.N += 1
            self.plot_frame_set()
            self.update_N_lbl()
    
    def N_subtract(self):
        
        if self.N == 1:
            pass
        else:
            self.N -= 1
            self.plot_frame_set()
            self.update_N_lbl()
    
    def plot_frame_set(self):
    
        self.axU.clear()
        self.axD.clear()
        
        data = self.frame_data[self.N-1]
        
        #data = {d.attrib['name']:
        #        np.reshape(
        #            np.fromstring(d.text, sep=';'),
        #            (int(d.attrib['x']),int(d.attrib['y'])))
        #        for d in self.frames[self.N-1] if d.tag=='Data'}
                
        self.axU.imshow(data['IntensityUp'].reshape((self.frame_x, self.frame_y)))
        self.axD.imshow(data['IntensityDown'].reshape((self.frame_x, self.frame_y)))
        
        self.canvas.draw()
    
    def process_frame_element(self, i, f):
    
        new_dict = {d.attrib['name']:
                    np.fromstring(d.text, sep=';').reshape(
                        (self.frame_x, self.frame_y)
                        )
                    for d in f if d.tag == 'Data'}
        
        progress = (i+1)/self.nframes*100
        self.progress_bar.setValue(int(progress))
        
        return new_dict
    
    def readXML(self):
        
        self.N = 1
        self.tree = ET.parse(self.filepath)
        self.root = self.tree.getroot()
        self.frames = [e for e in self.root if e.tag == 'Frame']
        self.nframes = len(self.frames)
        
        self.frame_x = int(self.frames[0].find('Data').attrib['x'])
        self.frame_y = int(self.frames[0].find('Data').attrib['y'])
        
        self.frame_data = [self.process_frame_element(i, f)
                           for i,f in enumerate(self.frames)]

        self.plot_frame_set()
        self.update_N_lbl()
        
    def start_frame_video(self):
        
        if self.N == self.nframes:
            pass
        else:
            self.video_thread = RunFrameVideo()
            self.video_thread.new_frame.connect(self.plot_next_image)
            self.video_thread.start()
            
    def plot_next_image(self):
    
        if self.N == self.nframes:
            self.stop_frame_video()
        else:
            self.N_add()
            self.plot_frame_set()

    def stop_frame_video(self):
    
        self.video_thread.run_video = False

#class PlayWorker(QObject):
#
#    finished = Signal()
#    ready = Signal()
#    
#    def __init__(self):
#    
#        self.run_video = True
#    
#    def new_frame(self):
#        
#        self.run_video = True
#        self.finished.emit()
        
    
class RunFrameVideo(QThread):
    
    new_frame = Signal()
    
    def __init__(self):
        
        QThread.__init__(self)
        
        self.run_video = True
        
    def __del__(self):
        
        self.wait()
    
    def run(self):
        
        while self.run_video:
            
            self.new_frame.emit()
            time.sleep(0.5)
    
if __name__ == '__main__':

    app = QApplication(sys.argv)
    w = SixT2Visual()
    sys.exit(app.exec_())