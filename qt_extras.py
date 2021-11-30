from PyQt5.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QPushButton
from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

class PlottingWindow(QWidget):

    def __init__(self):
    
        super(PlottingWindow, self).__init__()
        
        self.layout = QVBoxLayout()
        
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.tools = NavigationToolbar(self.canvas, self)
        self.ax = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
        
        self.layout.addWidget(self.canvas)
        
        self.tool_lo = QHBoxLayout()
        self.tool_lo.addWidget(self.tools)
        self.tool_lo.addStretch()
        
        self.reset_axes_btn = QPushButton('Reset axes')
        self.reset_axes_btn.clicked.connect(self.reset_axes)
        self.tool_lo.addWidget(self.reset_axes_btn)
        self.layout.addLayout(self.tool_lo)
        
        self.setLayout(self.layout)
        
    def reset_axes(self):
       
       s = 0
       if len(self.ax.lines)<1: pass
       else:
           while True:
               start = self.ax.lines[s]
               if len(start.get_xdata())<1:
                   s += 1
               else:
                   break
           
           x = start.get_xdata()
           y = start.get_ydata()
           
           new_x = [x.min(), x.max()]
           new_y = [y.min(), y.max()]
           
           for i in range(s, len(self.ax.lines)):
               x = self.ax.lines[i].get_xdata()
               y = self.ax.lines[i].get_ydata()
               
               if len(x)>1 and len(y)>1:
                   if x.min()<new_x[0]: new_x[0] = x.min()
                   if x.max()>new_x[1]: new_x[1] = x.max()
                   if y.min()<new_y[0]: new_y[0] = y.min()
                   if y.max()>new_y[1]: new_y[1] = y.max()
           
           if new_x[0] == new_x[1]:
               new_x[0] -= 0.5
               new_x[1] += 0.5
           if new_y[0] == new_y[1]:
               new_y[0] -= 0.5
               new_y[1] += 0.5
               
           self.ax.set_xlim(new_x[0]-0.05*(new_x[1]-new_x[0]),new_x[1]+0.05*(new_x[1]-new_x[0]))
           self.ax.set_ylim(new_y[0]-0.05*(new_y[1]-new_y[0]),new_y[1]+0.05*(new_y[1]-new_y[0]))
           
           self.canvas.draw()