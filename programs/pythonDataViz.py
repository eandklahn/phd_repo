import xml.etree.ElementTree as ET
import sys
import xml
import time

# EXTERNAL
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
from matplotlib.ticker import MaxNLocator

from PyQt5.QtGui import QIcon, QFont
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QWidget, QApplication, QPushButton, QGridLayout, QLabel, QComboBox, QStackedWidget, QSlider,
                             QDoubleSpinBox, QFormLayout, QCheckBox, QSpinBox, QVBoxLayout, QHBoxLayout, QMessageBox, QFileDialog)
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    
class DataViz(QWidget):


    def __init__(self):
        
        super().__init__()
        
        self.initUI()
    
    def initUI(self):
        
        # Opening a QWidget with a box layout
        LO = QVBoxLayout()
        
        # Adding the layout for file-related stuff
        fileLO = QHBoxLayout()
        
        openFilebtn = QPushButton('&Open')
        openFilebtn.clicked.connect(self.openFile)
        fileLO.addWidget(openFilebtn)
        
        self.currentFilelbl = QLabel()
        fileLO.addWidget(self.currentFilelbl)
        
        fileLO.addStretch()
        LO.addLayout(fileLO)
        
        # Adding the layout for frame-related stuff
        frameLO = QHBoxLayout()
        
        self.showFrameInfobtn = QPushButton('&Frame info')
        self.showFrameInfobtn.clicked.connect(self.displayFrameInfo)
        frameLO.addWidget(self.showFrameInfobtn)
        
        self.frameNumberBox = QSpinBox()
        self.frameNumberBox.setMinimum(1)
        self.frameNumberBox.valueChanged.connect(self.plotImage)
        frameLO.addWidget(self.frameNumberBox)
        
        self.noOfFrameslbl = QLabel()
        frameLO.addWidget(self.noOfFrameslbl)
        
        frameLO.addStretch()
        LO.addLayout(frameLO)
        
        # Adding the layout for image-related stuff
        
        sliderLO = QHBoxLayout()
        self.minvalueSlider = QSlider(Qt.Horizontal)
        self.minvalueSlider.sliderMoved.connect(self.plotImage)
        sliderLO.addWidget(self.minvalueSlider)
        
        self.maxvalueSlider = QSlider(Qt.Horizontal)
        self.maxvalueSlider.sliderMoved.connect(self.plotImage)
        sliderLO.addWidget(self.maxvalueSlider)
        
        self.areaDetectorImages = Figure()
        self.dynamic_canvas = FigureCanvas(self.areaDetectorImages)
        self.toolbar = NavigationToolbar(self.dynamic_canvas, self)
        
        self.axU = self.areaDetectorImages.add_subplot(121)
        self.axD = self.areaDetectorImages.add_subplot(122)
        
        # Containers for the objects that are returned by the imshow command
        self.imageUobject = None
        self.imageDobject = None
        
        LO.addWidget(self.dynamic_canvas)
        
        # Adding controls for the plotting
        controlsLO = QHBoxLayout()
        self.colormapBox = QComboBox()
        self.colormapBox.addItems(['hot', 'viridis', 'inferno', 'plasma', 'Greys'])
        self.colormapBox.currentIndexChanged.connect(self.plotImage)
        
        controlsLO.addWidget(self.toolbar)
        controlsLO.addWidget(self.colormapBox)
        
        controlsLO.addStretch()
        LO.addLayout(controlsLO)
        
        LO.addLayout(sliderLO)
        LO.addStretch()
        
        self.path = None
        self.tree = None
        self.root = None
        self.xmlInfo = None
        self.dataUArray = None
        self.dataDArray = None
        
        # Setting layout and opening the widget
        self.setLayout(LO)
        self.show()
    
    def displayFrameInfo(self):
        
        if self.root is not None:
            
            frameNo = self.frameNumberBox.value()
            idx = self.xmlInfo['framestartIndex']+frameNo-1
            
            msg = QMessageBox()
            msg.setStandardButtons(QMessageBox.Ok)
            msg.setWindowTitle('Frame info for #{}'.format(frameNo))
            
            frame = self.root[idx]
            
            infostr = ''
            info = frame.attrib
            for elem in frame.iter():
                if elem.tag != 'Data':
                    info.update(elem.attrib)
                    info.update({elem.tag: elem.text})
            
            for key, value in info.items():
                infostr += '{}: {}\n'.format(key, value)
            
            msg.setText(infostr)
            msg.exec_()
            
            
            
    def openFile(self):
        
        self.path = QFileDialog.getOpenFileName()[0]
        
        try:
            
            self.tree = ET.parse(self.path)
            self.root = self.tree.getroot()

            self.xmlInfo = self.readXMLinfo(self.root)
            
            self.frameNumberBox.setMaximum(len(self.root)-self.xmlInfo['framestartIndex'])
            self.maxvalueSlider.setValue(self.maxvalueSlider.maximum())
            
            self.currentFilelbl.setText(self.path)
            self.noOfFrameslbl.setText('/{}'.format(len(self.root)-self.xmlInfo['framestartIndex']))
            
        except (FileNotFoundError, xml.etree.ElementTree.ParseError) as e:

            msg = QMessageBox()
            msg.setText('Something happened! Try again')
            msg.setStandardButtons(QMessageBox.Ok)
            details = ''
            if isinstance(e, xml.etree.ElementTree.ParseError): details = 'Selected file could not be parsed as an XML-file'
            elif isinstance(e, FileNotFoundError): details = 'The selected file could not be found'
            msg.setDetailedText(details)
            msg.exec_()
        
        finally:
            self.frameNumberBox.setValue(1)
            self.plotImage()
            
    def plotImage(self):
        
        self.axU.clear()
        self.axD.clear()
        
        try:
            frameNo = self.frameNumberBox.value()
            
            dataU = self.root[self.xmlInfo['framestartIndex']+frameNo-1][self.xmlInfo['datastartIndex']+0].text
            dataUArray = np.fromstring(dataU, sep=';')
            self.dataUArray = dataUArray.reshape((self.xmlInfo['y'],self.xmlInfo['x']))
            
            dataD = self.root[self.xmlInfo['framestartIndex']+frameNo-1][self.xmlInfo['datastartIndex']+1].text
            dataDArray = np.fromstring(dataD, sep=';')
            self.dataDArray = dataDArray.reshape((self.xmlInfo['y'],self.xmlInfo['x']))
            
            self.maxvalueSlider.setMaximum(max(np.amax(self.dataDArray), np.amax(self.dataUArray)))
            self.minvalueSlider.setMaximum(self.maxvalueSlider.value())
            
            
            self.imageUobject = self.axU.imshow(self.dataUArray, extent=[0,self.xmlInfo['x'],0,self.xmlInfo['y']],
                                                            aspect=self.xmlInfo['x']/self.xmlInfo['y'],
                                                            vmin=self.minvalueSlider.value(),
                                                            vmax=self.maxvalueSlider.value(),
                                                            cmap=self.colormapBox.currentText())
            
            self.imageDobject = self.axD.imshow(self.dataDArray, extent=[0,self.xmlInfo['x'],0,self.xmlInfo['y']],
                                                            aspect=self.xmlInfo['x']/self.xmlInfo['y'],
                                                            vmin=self.minvalueSlider.value(),
                                                            vmax=self.maxvalueSlider.value(),
                                                            cmap=self.colormapBox.currentText())
            
            self.dynamic_canvas.draw()
            
        except TypeError:
            None
            
        
    def readXMLinfo(self, root):

        xmlInfo = {}
        xmlInfo['framestartIndex'] = 0
        while root[xmlInfo['framestartIndex']].tag != 'Frame':
            xmlInfo['framestartIndex'] += 1
        
        xmlInfo['datastartIndex'] = 0
        while root[xmlInfo['framestartIndex']][xmlInfo['datastartIndex']].tag != 'Data':
            xmlInfo['datastartIndex'] += 1
            
        xmlInfo['x'] = int(root[xmlInfo['framestartIndex']][xmlInfo['datastartIndex']].attrib['x'])
        xmlInfo['y'] = int(root[xmlInfo['framestartIndex']][xmlInfo['datastartIndex']].attrib['y'])
        
        return xmlInfo
    
    def get_minvalueSlider(self):
        
        v = self.minvalueSlider.value()
        
        return v
    
    def get_maxvalueSlider(self):
        
        v = self.maxvalueSlider.value()
        
        return v
    
    def keyPressEvent(self, event):
        
        if event.key() == QtCore.Qt.Key_O:
            N = self.frameNumberBox.value()
            self.frameNumberBox.setValue(N-1)
        elif event.key() == QtCore.Qt.Key_P:
            N = self.frameNumberBox.value()
            self.frameNumberBox.setValue(N+1)
            event.accept()
            
    
if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = DataViz()
    sys.exit(app.exec_())