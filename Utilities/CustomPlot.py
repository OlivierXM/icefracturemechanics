"""
    Custom Plotting Functionality from Matplotlib

    Provides a Matlab-esque interface
    Tasks:
    - Add support for arbitrary **kwargs for most functions

    From https://github.com/OlivierXM/dolfinxTools
"""
import matplotlib.pyplot as plt
from matplotlib import font_manager as fontP
import math

## Inherited matplotlib.font_manager
# Create a custom font set
# - **kwargs (size = 12, weight = "bold", ....) 
class Font(fontP.FontProperties):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

## CustomPlot.PrebuiltPlot()
# Provide MatLab-esque plotting args for arbitrary set of axes
class PrebuiltPlot():
    ## PrebuiltPlot.__init__()
    # Instantiate the figure with subplots
    # - a, b: The row and columns in figure
    # - fontProp: The CustomPlot.Font set to use unless otherwise stated
    def __init__(self, a, b, fontProp, shareX = True, **kwargs):
        self._font = fontProp
        self._fig, (self._axes) = plt.subplots(a, b, sharex = shareX, **kwargs)
        if (a*b == 1):
            self._axes = [self._axes, self._axes]

        self._row = a
        self._col = b
        self._size = a*b    

    ## PrebuiltPlot.plot()
    # Plot data on the specified axis
    # - axisN: The zero-based axis index
    # - datax, datay: np.arrays to plot
    # - **kwargs
    def plot(self, axisN, datax, datay, styleIn = '-', labelIn = None, lwIn = 2):
        self.convertN(axisN).plot(datax, datay, styleIn, label = labelIn, lw = lwIn)

    def xlabel(self, axisN, strIn, fontIn = None):
        if (fontIn == None): fontIn = self._font

        if (axisN < 0):
            for i in range(self._size):
                self.convertN(i).set_xlabel(strIn, fontproperties = fontIn)
        else:
            self.convertN(axisN).set_xlabel(strIn, fontproperties = fontIn)
   
    def ylabel(self, axisN, strIn, fontIn = None):
        if (fontIn == None): fontIn = self._font
            
        if (axisN < 0):
            for i in range(self._size):
                self.convertN(i).set_ylabel(strIn, fontproperties = fontIn)
        else:
            self.convertN(axisN).set_ylabel(strIn, fontproperties = fontIn)

    def legend(self, axisN, locIn = 'lower center', fontIn = None, **kwargs):
        if not(isinstance(axisN, list)): axisN = [axisN]
        if (fontIn == None): fontIn = self._font
        for axisI in axisN:
            self.convertN(axisI).legend(loc = locIn, prop = fontIn, **kwargs)

    def axis(self, axisN, xl=None, xu=None, yl=None, yu=None):
        if (axisN < 0):
            for i in range(self._size):
                self.convertN(i).set(xlim=(xl, xu), ylim=(yl, yu))
        else:
            self.convertN(axisN).set(xlim=(xl, xu), ylim=(yl, yu))

    def grid(self, axisN, whichIn = 'major'):
        if not(isinstance(axisN, list)):
            axisN = [axisN]

        for axisI in axisN:
            self.convertN(axisI).grid(which = whichIn)
            

    def title(self, strIn, fontIn = None):
        if (fontIn == None): fontIn = self._font
        self._fig.suptitle(strIn, fontproperties = fontIn)

    def show(self, blockIn = False):
        plt.show(block = blockIn)
        waitIt = input("Press enter to continue:.....")
        


    ## PrebuiltPlot.convertN()
    # Return the zero-based nxn or nx1 or 1xn axis i.e.
    # 2x2 => [0 2]
    #        [1 3]
    # - N: The zero-based index of the axis
    def convertN(self, N):
        floorC = math.floor(N / self._row)
        floorR = (N - floorC * self._row)
        if (self._row == 1 or self._col== 1):
            return self._axes[N]
        else:
            return self._axes[floorR, floorC]