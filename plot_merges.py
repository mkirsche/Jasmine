'''
A simple visualization for comparing Jasmine results to those of SURVIVOR on human chr1
It assumes that the exact lines and points to plot, along with their colors, have already
been determined with the companion program src/VisualizationPrep.java.

This program takes in a single command line argument - the name of the file with points/lines to plot
'''

# Lots of matplotlib imports - we need Qt5Agg for scrollbar
import matplotlib
matplotlib.use('Qt5Agg')

# Matplotlib's libraries
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import collections as mc
from matplotlib.lines import Line2D

# Other imports
import numpy as np
import pylab as pl
import sys
from PyQt5 import QtWidgets, QtCore

# No yellow because it's hard to see so use brown instead
colors = ['red', 'orange', 'brown', 'green', 'blue', 'purple', 'pink', 'gray', 'black']

xs = [] # x values are genomic positions of variants
ys = [] # y values are the sample ID
cs = [] # Color of each point based on its variant type
xline = [] # x-coordinate pairs for merged variants
yline = [] # y-coordinate pairs for merged variants
linecs = [] # Color of each line segment based on which software merged that pair


'''
Each line will contain either a point (x, y) or a line segment (x1, y1, x2, y2).
Both types of lines have the option of adding a color value afterwards.
'''
with open(sys.argv[1], "r") as f:
  for line in f.readlines():
    tokens = line.split()
    if len(tokens) == 2: # Point with no color
        xs.append(int(tokens[0]))
        ys.append(int(tokens[1]))
        cs.append(colors[ys[len(ys)-1]])
    elif len(tokens) == 4: # Line segment with no color
        xline.append([int(tokens[0]), int(tokens[2])])
        yline.append([int(tokens[1]), int(tokens[3])])
        linecs.append('black')
    elif len(tokens) == 3: # Point with color
        xs.append(int(tokens[0]))
        ys.append(int(tokens[1]))
        cs.append(colors[int(tokens[2])])
    elif len(tokens) == 5: # Line segment with color
        xline.append([int(tokens[0]), int(tokens[2])])
        yline.append([int(tokens[1]), int(tokens[3])])
        linecs.append(colors[2*int(tokens[4])]) # Double color value so not too similar

#plt.scatter(xs, ys, c = cs)
#for i in range(0, len(xline)):
#  plt.plot(xline[i], yline[i], c = linecs[i])

# A window to show a plot with scrolling along the x-axis enabled
class ScrollableWindow(QtWidgets.QMainWindow):

    # Here step is what proportion of x-axis to show at once
    def __init__(self, fig, ax, step=0.01):
        plt.close("all")
        if not QtWidgets.QApplication.instance():
            self.app = QtWidgets.QApplication(sys.argv)
        else:
            self.app = QtWidgets.QApplication.instance() 

        QtWidgets.QMainWindow.__init__(self)
        self.widget = QtWidgets.QWidget()
        self.setCentralWidget(self.widget)
        self.widget.setLayout(QtWidgets.QVBoxLayout())
        self.widget.layout().setContentsMargins(0,0,0,0)
        self.widget.layout().setSpacing(0)

        self.fig = fig
        self.ax = ax
        self.canvas = FigureCanvas(self.fig)
        self.canvas.draw()
        self.scroll = QtWidgets.QScrollBar(QtCore.Qt.Horizontal)
        self.step = step
        self.setupSlider()
        self.nav = NavigationToolbar(self.canvas, self.widget)
        self.widget.layout().addWidget(self.nav)
        self.widget.layout().addWidget(self.canvas)
        self.widget.layout().addWidget(self.scroll)

        self.canvas.draw()
        self.show()
        self.app.exec_()

    def setupSlider(self):
        self.lims = np.array(self.ax.get_xlim())
        self.scroll.setPageStep(self.step*100)
        self.scroll.actionTriggered.connect(self.update)
        self.update()

    # Update the window limits based on the scrollbar position
    def update(self, evt=None):
        r = self.scroll.value()/((1+self.step)*100)
        l1 = self.lims[0]+r*np.diff(self.lims)
        l2 = l1 +  np.diff(self.lims)*self.step
        self.ax.set_xlim(l1,l2)
        self.fig.canvas.draw_idle()

fig, ax = plt.subplots()

# Set the x-axis to go from 0 to the last variant position
plt.xlim(0, max(xs))

# Plot the variant points
ax.scatter(xs, ys, c = cs)

# Add axis labels and title
ax.set_ylabel('Sample ID')
ax.set_xlabel('Position (chr1)')
ax.set_yticks(np.arange(0, max(ys)+1))
ax.set_title('chr1')

# Plot the line segments
for i in range(0, len(xline)):
  if linecs[i] == colors[4]:
    ax.plot(xline[i], yline[i], c = linecs[i], linestyle='dotted')
  else:
    ax.plot(xline[i], yline[i], c = linecs[i])
custom_lines = [Line2D([0], [0], color=colors[2], lw=4),
                Line2D([0], [0], color=colors[4], lw=4),
                Line2D([0], [0], color=colors[6], lw=4)]

# Add legend for merging software colors
legend1 = plt.legend(custom_lines, ['Jasmine', 'SURVIVOR', 'BOTH'], bbox_to_anchor=(.3, 1.05), ncol = 3)
ax.add_artist(legend1)

# Add legend for variant type colors
patches = [mpatches.Patch(color=colors[0], label='INS'), mpatches.Patch(color=colors[1], label='DEL'),
    mpatches.Patch(color=colors[2], label='DUP'), mpatches.Patch(color=colors[3], label='INV')]
legend2 = plt.legend(handles=patches, bbox_to_anchor=(.9, 1.05), ncol=len(patches))
ax.add_artist(legend2)

# Generate the plot with a scrolling window
a = ScrollableWindow(fig,ax)


