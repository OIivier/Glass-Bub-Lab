
#this is a script designed to streamline the analysis of redshirt imaging data files. There is a GUI that allows users to select data files, read them, identify areas of high activity, browse through time-series create IBI plots and return maps and create figures in several file formats.
#Written by Alex Hodge
#last updated Oct 28, 2008

import matplotlib
matplotlib.use('WX')
from matplotlib.backends.backend_wx import Toolbar, FigureCanvasWx,\
     FigureManager

import numpy
from numpy import arange, sin, pi

from matplotlib.figure import Figure
from matplotlib.axes import Subplot

from wxPython.wx import *

import struct

from pylab import *
from time import *

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg 
from matplotlib.backends.backend_wx import NavigationToolbar2Wx

from matplotlib.backends.backend_wx import _load_bitmap
from matplotlib.figure import Figure
#from matplotlib.numerix.mlab import rand
#from matplotlib.numerix.mlab import mean

import math
from scipy import signal, arange, sin, pi, linspace, transpose
from scipy.signal import buttord, butter, lfilter
import random

import wx
import os,sys


coords=[]

global f
global f1
global FileName
FileString=[]

print '\n'
print 'This program reads large data files from RedShirt cameras and allows \nfor pixel selection for further processing\n'

HeaderLength=2560;
snapshot = []
snapvec=[]
snapvecReorg=[];
rangevec = []
snapshotReorg = []
rangevecReorg = []
pixelvec=[]
LongTSPixelMat=[]
FilteredLongTS=[]
DataTime=localtime()
Threshold=[0,0, 0, 0, 0]
ibis=[]
TSTemporaryVec=sys.argv[0]
ChosenIBIMap=0
RetMapLower=1
RetMapUpper=10000
RetMapLims=[0,0]
ChosenIBIMap=[]



QSize=40 #used so I can dynamically set the quadrant size
print 'Creating map vector'
mapvec=[]
for i in range(QSize**2):
        app=4*i+1
        mapvec.append(app)
        if (i+1)%QSize==0:
                for j in range(QSize):
                        app=4*QSize-2-4*j+4*QSize*(i/QSize)
                        mapvec.append(app)

for i in range(QSize**2):
        app=4*i+3  +  4*QSize*(QSize-2*(i/QSize)-1)
        mapvec.append(app)
        if (i+1)%QSize==0:
                for j in range(QSize):
                        app=4*QSize-4*j+4*QSize*(i/QSize)+4*QSize*(QSize-2*(i/QSize)-1)
                        mapvec.append(app)
print 'done\n'

exx=68
wyy=39
ind=(80*(79-wyy)+exx)
print mapvec[ind]-1
exx=1
wyy=79
ind=80*(79-wyy)+exx
print mapvec[ind]
print mapvec[1]
threshold=DataTime


print 'Initialising data matrix...'
dummy=[[]]
for i in range(6399):
	dummy.append([])
print 'done\n'

if(threshold[1]>10):
	os.system("rm "+TSTemporaryVec)

class MyDialoga(wx.Dialog):# Not used now
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, size=(350,300))

        sizer =  self.CreateTextSizer('My Buttons')
        sizer.Add(wx.Button(self, -1, 'Button'), 0, wx.ALL, 5)
        sizer.Add(wx.Button(self, -1, 'Button'), 0, wx.ALL, 5)
        sizer.Add(wx.Button(self, -1, 'Button'), 0, wx.ALL, 5)
        sizer.Add(wx.Button(self, -1, 'Button'), 0, wx.ALL|wx.ALIGN_CENTER, 5)
        sizer.Add(wx.Button(self, -1, 'Button'), 0, wx.ALL|wx.EXPAND, 5)
        self.SetSizer(sizer)

class MyNavigationToolbar(NavigationToolbar2WxAgg):
    """
    Extend the default wx toolbar with your own event handlers
    """
    ON_CUSTOM = wxNewId()
    def __init__(self, canvas, cankill):
        NavigationToolbar2WxAgg.__init__(self, canvas)
        self.AddSimpleTool(self.ON_CUSTOM, _load_bitmap('stock_left.xpm'),
                           'Click me', 'Activate custom contol')
        EVT_TOOL(self, self.ON_CUSTOM, self._on_custom)

    def _on_custom(self, evt):
        ax = self.canvas.figure.axes[0]

        # generate a random location can color
        x,y = tuple(rand(2))
        rgb = tuple(rand(3))

        # add the text and draw
        ax.text(x, y, 'You clicked me',
                transform=ax.transAxes,
                color=rgb)
        self.canvas.draw()
        evt.Skip()

#####################################

class CanvasFrameMap(wx.Frame):

    def __init__(self, ):
        wx.Frame.__init__(self,None,-1,
                         'Activity Map',size=(750,750))

        self.SetBackgroundColour(wx.NamedColor("BLACK"))

        self.figure = Figure(figsize=(5,5))
        self.axes = self.figure.add_subplot(111)


	duder=array(rangevecReorg)
	duder.shape = 80,80


        self.axes.imshow(duder, origin='upper', cmap=cm.hot, extent=(0,80,0,80))
        self.axes.set_xlabel(str(FileString[0]))
        self.axes.set_ylabel('y')
        self.figure_canvas = FigureCanvas(self, -1, self.figure)

        # Note that event is a MplEvent
        self.figure_canvas.mpl_connect('motion_notify_event', self.UpdateStatusBar)
        self.figure_canvas.Bind(wx.EVT_ENTER_WINDOW, self.ChangeCursor)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.figure_canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

        self.statusBar = wx.StatusBar(self, -1)
        self.statusBar.SetFieldsCount(1)
        self.SetStatusBar(self.statusBar)

        self.toolbar = NavigationToolbar2Wx(self.figure_canvas)
        self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.toolbar.Show()

    def ChangeCursor(self, event):
        self.figure_canvas.SetCursor(wx.StockCursor(wx.CURSOR_BULLSEYE))

    def UpdateStatusBar(self, event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            self.statusBar.SetStatusText(( "x= " + str(int(x)) +
                                           "  y=" +str(int(y)) ),
                                           0)


########################################


class CanvasFrameTS(wxFrame):
    
    def __init__(self):
        wxFrame.__init__(self,None,-1,
                         'Time Series',size=(1450,1100))

        self.SetBackgroundColour(wxNamedColor("BLACK"))



        self.figure = Figure(figsize=(5,5), dpi=100)
        self.axes1 = self.figure.add_subplot(511)
	self.axes2 = self.figure.add_subplot(512)
        self.axes3 = self.figure.add_subplot(513)
	self.axes4 = self.figure.add_subplot(514)
        self.axes5 = self.figure.add_subplot(515)
        self.axes5.set_xlabel(str(FileString[0]), fontsize = 8)     
   
	if (len(pixelvec)>0):###Might try all on one axis system
	        self.axes1.plot(FilteredLongTS[0])
		tvec=Threshold[0]+ones(len(FilteredLongTS[0]))
		self.axes1.plot(tvec)

	if (len(pixelvec)>1):
        	self.axes2.plot(FilteredLongTS[1])
		tvec=Threshold[1]+ones(len(FilteredLongTS[0]))
		self.axes2.plot(tvec)

	if (len(pixelvec)>2):
	        self.axes3.plot(FilteredLongTS[2])
		tvec=Threshold[2]+ones(len(FilteredLongTS[0]))
		self.axes3.plot(tvec)

	if (len(pixelvec)>3):
        	self.axes4.plot(FilteredLongTS[3])
		tvec=Threshold[3]+ones(len(FilteredLongTS[0]))
		self.axes4.plot(tvec)

	if (len(pixelvec)>4):
        	self.axes5.plot(FilteredLongTS[4])
		tvec=Threshold[4]+ones(len(FilteredLongTS[0]))
		self.axes5.plot(tvec)

        self.canvas = FigureCanvas(self, -1, self.figure)

        self.sizer = wxBoxSizer(wxVERTICAL)
        self.sizer.Add(self.canvas, 1, wxTOP | wxLEFT | wxEXPAND)
        # Capture the paint message        
        EVT_PAINT(self, self.OnPaint)        

        self.toolbar = MyNavigationToolbar(self.canvas, True)
        self.toolbar.Realize()
        if wxPlatform == '__WXMAC__':
            # Mac platform (OSX 10.3, MacPython) does not seem to cope with
            # having a toolbar in a sizer. This work-around gets the buttons
            # back, but at the expense of having the toolbar at the top
            self.SetToolBar(self.toolbar)
        else:
            # On Windows platform, default window size is incorrect, so set
            # toolbar width to figure width.
            tw, th = self.toolbar.GetSizeTuple()
            fw, fh = self.canvas.GetSizeTuple()
            # By adding toolbar in sizer, we are able to put it at the bottom
            # of the frame - so appearance is closer to GTK version.
            # As noted above, doesn't work for Mac.
            self.toolbar.SetSize(wxSize(fw, th))
            self.sizer.Add(self.toolbar, 0, wxLEFT | wxEXPAND)

        # update the axes menu on the toolbar
        self.toolbar.update()  
        self.SetSizer(self.sizer)
        self.Fit()


    def OnPaint(self, event):
        self.canvas.draw()
        event.Skip()

class CanvasFrameIBIs(wxFrame):
    
    def __init__(self):
        wxFrame.__init__(self,None,-1,
                         'IBIs',size=(1450,1100))


	IBIs=[]
	IBITimes=[]
	for i in range(len(pixelvec)):

		ibis=[]
		timeZ=[]
		j=0
		thresh=Threshold[i]
		print thresh
		xlp=FilteredLongTS[i]
		for i in range(2*71987):
		#for i in range(len(longvec)-1):
			j+=1
			if ((xlp[i]<thresh)&(xlp[i+1]>thresh)):
				ibis.append(j)
				timeZ.append(i)
				j=0
		print 'ibis length, i ',len(ibis),i
		IBIs.append(ibis)
		IBITimes.append(timeZ)


	#	print LongTSPixelMat
	print len(LongTSPixelMat)	
	

	#	print LongTSPixelMat
	print len(LongTSPixelMat)	




        self.SetBackgroundColour(wxNamedColor("BLACK"))

	print Threshold

        self.figure = Figure(figsize=(5,5), dpi=100)
        self.axes1 = self.figure.add_subplot(511)

	self.axes2 = self.figure.add_subplot(512)

        self.axes3 = self.figure.add_subplot(513)

	self.axes4 = self.figure.add_subplot(514)

        self.axes5 = self.figure.add_subplot(515)





        self.axes5.set_xlabel(str(FileString[0]), fontsize = 8)
	if (len(pixelvec)>0):
	        self.axes1.plot(IBITimes[0],IBIs[0],'.')
		self.axes1.text(10, 10, str(coords[0]), color='black', size=10)
        	self.axes1.set_xlim([-10, 144000])
        	self.axes1.set_ylim([-10, 300])
       	 	self.axes1.set_ylabel(str(coords[0]), fontsize = 8)

	if (len(pixelvec)>1):
        	self.axes2.plot(IBITimes[1],IBIs[1],'.')
		self.axes2.text(10, 10, str(coords[1]), color='black', size=10)
 	        self.axes2.set_xlim([-10, 144000])
 	        self.axes2.set_ylim([-10, 300])
        	self.axes2.set_ylabel(str(coords[1]), fontsize = 8)

	if (len(pixelvec)>2):
	        self.axes3.plot(IBITimes[2],IBIs[2],'.')
		self.axes3.text(10, 10, str(coords[2]), color='black', size=10)
	        self.axes3.set_xlim([-10, 144000])
	        self.axes3.set_ylim([-10, 300])
	        self.axes3.set_ylabel(str(coords[2]), fontsize = 8)

	if (len(pixelvec)>3):
        	self.axes4.plot(IBITimes[3],IBIs[3],'.')
		self.axes4.text(10, 10, str(coords[3]), color='black', size=10)
	        self.axes4.set_xlim([-10, 144000])
	        self.axes4.set_ylim([-10, 300])
	        self.axes4.set_ylabel(str(coords[3]), fontsize = 8)
	
	if (len(pixelvec)>4):
        	self.axes5.plot(IBITimes[4],IBIs[4],'.')
		self.axes5.text(10, 10, str(coords[4]), color='black', size=10)
	        self.axes5.set_xlim([-10, 144000])
	        self.axes5.set_ylim([-10, 300])
	        self.axes5.set_ylabel(str(coords[4]), fontsize = 8)

        self.canvas = FigureCanvas(self, -1, self.figure)

        self.sizer = wxBoxSizer(wxVERTICAL)
        self.sizer.Add(self.canvas, 1, wxTOP | wxLEFT | wxEXPAND)
        # Capture the paint message        
        EVT_PAINT(self, self.OnPaint)        

        self.toolbar = MyNavigationToolbar(self.canvas, True)
        self.toolbar.Realize()
        if wxPlatform == '__WXMAC__':
            # Mac platform (OSX 10.3, MacPython) does not seem to cope with
            # having a toolbar in a sizer. This work-around gets the buttons
            # back, but at the expense of having the toolbar at the top
            self.SetToolBar(self.toolbar)
        else:
            # On Windows platform, default window size is incorrect, so set
            # toolbar width to figure width.
            tw, th = self.toolbar.GetSizeTuple()
            fw, fh = self.canvas.GetSizeTuple()
            # By adding toolbar in sizer, we are able to put it at the bottom
            # of the frame - so appearance is closer to GTK version.
            # As noted above, doesn't work for Mac.
            self.toolbar.SetSize(wxSize(fw, th))
            self.sizer.Add(self.toolbar, 0, wxLEFT | wxEXPAND)

        # update the axes menu on the toolbar
        self.toolbar.update()  
        self.SetSizer(self.sizer)
        self.Fit()


    def OnPaint(self, event):
        self.canvas.draw()
        event.Skip()

class CanvasFrameRetMap(wxFrame):
    
    def __init__(self):
        wxFrame.__init__(self,None,-1,
                         'Return Map',size=(1700,900))

	IBIs=[]
	IBITimes=[]
	for i in range(len(pixelvec)):

		ibis=[]
		timeZ=[]
		j=0
		thresh=Threshold[i]
		print thresh
		xlp=FilteredLongTS[i]

		for i in range(2*71990):
		#for i in range(len(longvec)-1):
			j+=1
			if ((xlp[i]<thresh)&(xlp[i+1]>thresh)):
				ibis.append(j)
				timeZ.append(i)
				j=0
		print 'ibis length, i ',len(ibis),i
		IBIs.append(ibis)
		IBITimes.append(timeZ)
	

	#	print LongTSPixelMat
	print len(LongTSPixelMat)	




	IBIsA=[]
	IBIsB=[]

	for i in range(len(pixelvec)):
		ibisa=IBIs[i][0:len(IBIs[i])-1]
		ibisb=IBIs[i][1:len(IBIs[i])]
		IBIsA.append(ibisa)
		IBIsB.append(ibisb)	

	print 'IBI len ',len(IBIs[0])

	print 'IBIsA len ',len(IBIsA)

        self.SetBackgroundColour(wxNamedColor("BLACK"))

#	print Threshold

        self.figure = Figure(figsize=(6,6), dpi=150)
        self.axes1 = self.figure.add_subplot(231)

	self.axes2 = self.figure.add_subplot(232)
	self.axes2.set_title(str(FileString[0]))

        self.axes3 = self.figure.add_subplot(233)

	self.axes4 = self.figure.add_subplot(234)

        self.axes5 = self.figure.add_subplot(235)





#        self.axes5.set_xlabel(str(FileString[0]), fontsize = 8)
        self.axes5.set_xlabel('dfgs', fontsize = 8)


	if (len(pixelvec)>0):
	        self.axes1.plot(IBIsA[0],IBIsB[0],'.',markersize=2)
		self.axes1.plot(range(300))
        	self.axes1.set_xlim([-10, 300])
        	self.axes1.set_ylim([-10, 300])
        	self.axes1.set_xlabel(str(coords[0]), fontsize = 8)

	if (len(pixelvec)>1):
	        self.axes2.plot(IBIsA[1],IBIsB[1],'.',markersize=2)
		self.axes2.plot(range(300))
        	self.axes2.set_xlim([-10, 300])
        	self.axes2.set_ylim([-10, 300])
        	self.axes2.set_xlabel(str(coords[1]), fontsize = 8)

	if (len(pixelvec)>2):
	        self.axes3.plot(IBIsA[2],IBIsB[2],'.',markersize=2)
		self.axes3.plot(range(300))
        	self.axes3.set_xlim([-10, 300])
        	self.axes3.set_ylim([-10, 300])
        	self.axes3.set_xlabel(str(coords[2]), fontsize = 8)

	if (len(pixelvec)>3):
	        self.axes4.plot(IBIsA[3],IBIsB[3],'.',markersize=2)
		self.axes4.plot(range(300))
        	self.axes4.set_xlim([-10, 300])
        	self.axes4.set_xlabel(str(coords[3]), fontsize = 8)
        	self.axes4.set_ylim([-10, 300])

	if (len(pixelvec)>4):
	        self.axes5.plot(IBIsA[4],IBIsB[4],'.',markersize=2)
		self.axes5.plot(range(300))
        	self.axes5.set_xlim([-10, 300])
        	self.axes5.set_ylim([-10, 300])
        	self.axes5.set_xlabel(str(coords[4]), fontsize = 8)

        self.canvas = FigureCanvas(self, -1, self.figure)

        self.sizer = wxBoxSizer(wxVERTICAL)
        self.sizer.Add(self.canvas, 1, wxTOP | wxLEFT | wxEXPAND)
        # Capture the paint message        
        EVT_PAINT(self, self.OnPaint)        

        self.toolbar = MyNavigationToolbar(self.canvas, True)
        self.toolbar.Realize()
        if wxPlatform == '__WXMAC__':
            # Mac platform (OSX 10.3, MacPython) does not seem to cope with
            # having a toolbar in a sizer. This work-around gets the buttons
            # back, but at the expense of having the toolbar at the top
            self.SetToolBar(self.toolbar)
        else:
            # On Windows platform, default window size is incorrect, so set
            # toolbar width to figure width.
            tw, th = self.toolbar.GetSizeTuple()
            fw, fh = self.canvas.GetSizeTuple()
            # By adding toolbar in sizer, we are able to put it at the bottom
            # of the frame - so appearance is closer to GTK version.
            # As noted above, doesn't work for Mac.
            self.toolbar.SetSize(wxSize(fw, th))
            self.sizer.Add(self.toolbar, 0, wxLEFT | wxEXPAND)

        # update the axes menu on the toolbar
        self.toolbar.update()  
        self.SetSizer(self.sizer)
        self.Fit()



    def OnPaint(self, event):
        self.canvas.draw()
        event.Skip()

class RetMapSection(wxFrame):
    
    
    def __init__(self):
        wxFrame.__init__(self,None,-1,
                         'Return Map',size=(1700,900))

	print RetMapLims[0]
	print RetMapLims[1]
	print ChosenIBIMap

	i = ChosenIBIMap[len(ChosenIBIMap)-1]-1
	k=i
	ibis=[]
	timeZ=[]
	


	j=0
	thresh=Threshold[i]
	print thresh
	xlp=FilteredLongTS[i]
	print len(xlp),'xlp'
	RMVec=xlp[RetMapLims[0]:RetMapLims[1]]
	print len(RMVec),'RMVec'
	for i in range(len(RMVec)-1):
	#for i in range(len(longvec)-1):
		j+=1
		if ((RMVec[i]<thresh)&(RMVec[i+1]>thresh)):
			ibis.append(j)
			timeZ.append(i)
			j=0
	print 'ibis length, i ',len(ibis),i



	ibisa=ibis[0:len(ibis)-1]
	ibisb=ibis[1:len(ibis)]


        self.SetBackgroundColour(wxNamedColor("BLACK"))

#	print Threshold

        self.figure = Figure(figsize=(6,6), dpi=150)
        self.axes1 = self.figure.add_subplot(111)
	self.axes1.set_title(str(FileString[0]),fontsize = 12)

        self.axes1.plot(ibisa,ibisb,'.',markersize=2)
	self.axes1.plot(range(300))
       	self.axes1.set_xlim([-10, 300])
       	self.axes1.set_ylim([-10, 300])


	print k
	print coords[k]
	print RetMapLims[0]
	print RetMapLims[1]
	

       	self.axes1.set_xlabel(str(coords[k])+'Frames '+str(RetMapLims[0])+'to '+str(RetMapLims[1]), fontsize = 8)

        self.canvas = FigureCanvas(self, -1, self.figure)

        self.sizer = wxBoxSizer(wxVERTICAL)
        self.sizer.Add(self.canvas, 1, wxTOP | wxLEFT | wxEXPAND)
        # Capture the paint message        
        EVT_PAINT(self, self.OnPaint)        

        self.toolbar = MyNavigationToolbar(self.canvas, True)
        self.toolbar.Realize()
        if wxPlatform == '__WXMAC__':
            # Mac platform (OSX 10.3, MacPython) does not seem to cope with
            # having a toolbar in a sizer. This work-around gets the buttons
            # back, but at the expense of having the toolbar at the top
            self.SetToolBar(self.toolbar)
        else:
            # On Windows platform, default window size is incorrect, so set
            # toolbar width to figure width.
            tw, th = self.toolbar.GetSizeTuple()
            fw, fh = self.canvas.GetSizeTuple()
            # By adding toolbar in sizer, we are able to put it at the bottom
            # of the frame - so appearance is closer to GTK version.
            # As noted above, doesn't work for Mac.
            self.toolbar.SetSize(wxSize(fw, th))
            self.sizer.Add(self.toolbar, 0, wxLEFT | wxEXPAND)

        # update the axes menu on the toolbar
        self.toolbar.update()  
        self.SetSizer(self.sizer)
        self.Fit()


    def OnPaint(self, event):
        self.canvas.draw()
        event.Skip()



	print RetMapLims[0]
	print RetMapLims[1]
	print ChosenIBIMap


class PlotFigure(wxFrame):
    def __init__(self):
        wxFrame.__init__(self, None, -1, "Test embedded wxFigure")

        self.fig = Figure((9,8), 75)
        self.canvas = FigureCanvasWx(self, -1, self.fig)
        self.toolbar = Toolbar(self.canvas)
        self.toolbar.Realize()

        # On Windows, default frame size behaviour is incorrect
        # you don't need this under Linux
        tw, th = self.toolbar.GetSizeTuple()
        fw, fh = self.canvas.GetSizeTuple()
        self.toolbar.SetSize(wxSize(fw, th))

        # Create a figure manager to manage things
        self.figmgr = FigureManager(self.canvas, 1, self)
        # Now put all into a sizer
        sizer = wxBoxSizer(wxVERTICAL)
        # This way of adding to sizer allows resizing
        sizer.Add(self.canvas, 1, wxLEFT|wxTOP|wxGROW)
        # Best to allow the toolbar to resize!
        sizer.Add(self.toolbar, 0, wxGROW)
        self.SetSizer(sizer)
        self.Fit()

    def plot_data(self):
        # Use ths line if using a toolbar
        a = self.fig.add_subplot(111)
        
        # Or this one if there is no toolbar
        #a = Subplot(self.fig, 111)
        
        t = numpy.arange(0.0,3.0,0.01)
        s = numpy.sin(2*numpy.pi*t)
        c = numpy.cos(2*numpy.pi*t)
        a.plot(t,s)
        a.plot(t,c)
        self.toolbar.update()

    def GetToolBar(self):
        # You will need to override GetToolBar if you are using an 
        # unmanaged toolbar in your frame
        return self.toolbar


class Slider(wx.Frame):
    def __init__(self, parent, id, title):

        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, (300, 150))
        panel = wx.Panel(self, -1)

        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.sld1 = wx.Slider(panel, -1, Threshold[0], -500, 500, wx.DefaultPosition, (250, -1),
                              wx.SL_AUTOTICKS | wx.SL_HORIZONTAL | wx.SL_LABELS)
        self.sld2 = wx.Slider(panel, -1, Threshold[1], -500, 500, wx.DefaultPosition, (250, -1),
                              wx.SL_AUTOTICKS | wx.SL_HORIZONTAL | wx.SL_LABELS)
        self.sld3 = wx.Slider(panel, -1, Threshold[2], -500, 500, wx.DefaultPosition, (250, -1),
                              wx.SL_AUTOTICKS | wx.SL_HORIZONTAL | wx.SL_LABELS)
        self.sld4 = wx.Slider(panel, -1, Threshold[3], -500, 500, wx.DefaultPosition, (250, -1),
                              wx.SL_AUTOTICKS | wx.SL_HORIZONTAL | wx.SL_LABELS)
        self.sld5 = wx.Slider(panel, -1, Threshold[4], -500, 500, wx.DefaultPosition, (250, -1),
                              wx.SL_AUTOTICKS | wx.SL_HORIZONTAL | wx.SL_LABELS)
        btn1 = wx.Button(panel, 8, 'Adjust')
        btn2 = wx.Button(panel, 9, 'Close')

        wx.EVT_BUTTON(self, 8, self.OnAdjust)
        wx.EVT_BUTTON(self, 9, self.OnClose)
        vbox.Add(self.sld1, 1, wx.ALIGN_CENTRE)
        vbox.Add(self.sld2, 1, wx.ALIGN_CENTRE)
        vbox.Add(self.sld3, 1, wx.ALIGN_CENTRE)
        vbox.Add(self.sld4, 1, wx.ALIGN_CENTRE)
        vbox.Add(self.sld5, 1, wx.ALIGN_CENTRE)
        hbox.Add(btn1, 1, wx.RIGHT, 10)
        hbox.Add(btn2, 1)
        vbox.Add(hbox, 0, wx.ALIGN_CENTRE | wx.ALL, 20)
        panel.SetSizer(vbox)


        self.SetSize((300, 400))

    def OnAdjust(self, event):
        val1 = self.sld1.GetValue()
	Threshold[0]=val1
        val2 = self.sld2.GetValue()
	Threshold[1]=val2
        val3 = self.sld3.GetValue()
	Threshold[2]=val3
        val4 = self.sld4.GetValue()
	Threshold[3]=val4
        val5 = self.sld5.GetValue()
	Threshold[4]=val5

	print val1
	print val2	
	print val3
	print val4
	print val5
	print Threshold
    def OnClose(self, event):
        self.Close()


class MyFrame(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(750, 450))

	self.FileName=''

        hbox  = wx.BoxSizer(wx.HORIZONTAL)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        vbox3 = wx.GridSizer(2,2,0,0)
        vbox4 = wx.BoxSizer(wx.VERTICAL)


        panel = wx.Panel(self, -1)
        pnl1 = wx.Panel(self, -1, style=wx.SIMPLE_BORDER)
        pnl2 = wx.Panel(self, -1, style=wx.SIMPLE_BORDER)
	pnl3 = wx.Panel(self, -1, style=wx.SIMPLE_BORDER)


        self.rb1 = wx.RadioButton(pnl3, -1, 'IBI Map 1', (10, 10), style=wx.RB_GROUP)
        self.rb2 = wx.RadioButton(pnl3, -1, 'IBI Map 2', (10, 30))
        self.rb3 = wx.RadioButton(pnl3, -1, 'IBI Map 3', (10, 50))
        self.rb4 = wx.RadioButton(pnl3, -1, 'IBI Map 4', (10, 70))
        self.rb5 = wx.RadioButton(pnl3, -1, 'IBI Map 5', (10, 90))
        self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb1.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb2.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb3.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb4.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb5.GetId())
        self.SetVal(True)
        self.sc1=wx.SpinCtrl(pnl3, -1, 'Start Frame', (15, 120), (80, -1), min=1, max=150000)
        self.sc2=wx.SpinCtrl(pnl3, -1, 'Stop Frame', (15, 150), (80, -1), min=1, max=150000)
        wx.Button(pnl3, 82, 'Limits', (10, 200))

        self.Bind(wx.EVT_BUTTON, self.OnLimits, id=82)

      	self.CreateStatusBar()
      	menuBar = wx.MenuBar()
      	menu = wx.Menu()
      	menu.Append(99,  "&Message Dialog", "Shows a Message Dialog")
      	menu.Append(100, "&Color Dialog", "Shows a Color Dialog")
      	menu.Append(101, "&Open File", "Shows a File Dialog")
      	menu.Append(102, "&Page Setup Dialog", "Shows a Page Setup Dialog")
      	menu.Append(103, "&Font Dialog", "Shows a Font Dialog")
      	menu.Append(104, "&Directory Dialog", "Shows a Directory Dialog")
      	menu.Append(105, "&SingleChoice Dialog", "Shows a SingleChoice Dialog")
      	menu.Append(106, "&TextEntry Dialog", "Shows a TextEntry Dialog")
      	menuBar.Append(menu, "&Do Something")
      	self.SetMenuBar(menuBar)


        self.lc = wx.ListCtrl(self, -1, style=wx.LC_REPORT)
        self.lc.InsertColumn(0, 'X')
        self.lc.InsertColumn(1, 'Y')
        self.lc.SetColumnWidth(0, 80)
        self.lc.SetColumnWidth(1, 80)
        vbox1.Add(pnl1, 1, wx.EXPAND | wx.ALL, 3)
        vbox1.Add(pnl2, 1, wx.EXPAND | wx.ALL, 3)
        vbox2.Add(self.lc, 1, wx.EXPAND | wx.ALL, 3)
        self.tc1 = wx.TextCtrl(pnl1, -1)
        self.tc2 = wx.TextCtrl(pnl1, -1)
        vbox3.AddMany([ (wx.StaticText(pnl1, -1, 'X'),0, wx.ALIGN_CENTER),
                        (self.tc1, 0, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL),
                        (wx.StaticText(pnl1, -1, 'Y'),0, wx.ALIGN_CENTER_HORIZONTAL),
                        (self.tc2,0)])
        pnl1.SetSizer(vbox3)
        vbox4.Add(wx.Button(pnl2, 10, 'Add'),   0, wx.ALIGN_CENTER| wx.TOP, 10)
        vbox4.Add(wx.Button(pnl2, 11, 'Remove'), 0, wx.ALIGN_CENTER|wx.TOP, 10)
        vbox4.Add(wx.Button(pnl2, 12, 'Clear'), 0, wx.ALIGN_CENTER| wx.TOP, 10)
        vbox4.Add(wx.Button(pnl2, 13, 'Process'), 0, wx.ALIGN_CENTER| wx.TOP, 10)
        pnl2.SetSizer(vbox4)

	hbox.Add(panel, 1, wx.EXPAND)
        hbox.Add(vbox1, 1, wx.EXPAND)
        hbox.Add(vbox2, 1, wx.EXPAND)
	hbox.Add(pnl3, 1, wx.EXPAND)
        self.SetSizer(hbox)



        self.Bind (wx.EVT_BUTTON, self.OnAdd, id=10)
        self.Bind (wx.EVT_BUTTON, self.OnRemove, id=11)
        self.Bind (wx.EVT_BUTTON, self.OnClear, id=12)
        self.Bind (wx.EVT_BUTTON, self.OnProcess, id=13)

      	self.Bind(wx.EVT_MENU, self.message, id=99)
      	self.Bind(wx.EVT_MENU, self.choosecolor, id=100)
      	self.Bind(wx.EVT_MENU, self.openfile, id=101)
      	self.Bind(wx.EVT_MENU, self.pagesetup, id=102)
      	self.Bind(wx.EVT_MENU, self.choosefont, id=103)
      	self.Bind(wx.EVT_MENU, self.opendir, id=104)
      	self.Bind(wx.EVT_MENU, self.singlechoice, id=105)
      	self.Bind(wx.EVT_MENU, self.textentry, id=106)

        wx.Button(panel, 1, 'Show IBIs', (40,50))
        self.Bind (wx.EVT_BUTTON, self.OnShowIBIs, id=1)
        wx.Button(panel, 2, 'Show Time Series', (40,100))
        self.Bind (wx.EVT_BUTTON, self.OnShowTS, id=2)
        wx.Button(panel, 3, 'Show Map', (40,150))
        self.Bind (wx.EVT_BUTTON, self.OnShowTSa, id=3)
        wx.Button(panel, 4, 'Adjust Thresholds', (40,200))
        self.Bind (wx.EVT_BUTTON, self.OnShowSlider, id=4)
        wx.Button(panel, 5, 'Read Data', (40,250))
        self.Bind (wx.EVT_BUTTON, self.OnReadData, id=5)
        wx.Button(panel, 6, 'Show RetMap', (40,300))
        self.Bind (wx.EVT_BUTTON, self.OnShowRetMap, id=6)


	self.SetStatusText("Current Data File: %s" % self.FileName)

    def OnLimits(self, event):
	RetMapLims[0]=self.sc1.GetValue()
	RetMapLims[1]=self.sc2.GetValue()
        frame = RetMapSection()
        frame.Show(true)
        frame.Centre()
	
    def SetVal(self, event):
        state1 = str(self.rb1.GetValue())
        state2 = str(self.rb2.GetValue())
        state3 = str(self.rb3.GetValue())
        state4 = str(self.rb4.GetValue())
        state5 = str(self.rb5.GetValue())

	if(self.rb1.GetValue()):
		ChosenIBIMap.append(1)
		print ChosenIBIMap
	if(self.rb2.GetValue()):
		ChosenIBIMap.append(2)
		print ChosenIBIMap
	if(self.rb3.GetValue()):
		ChosenIBIMap.append(3)
		print ChosenIBIMap
	if(self.rb4.GetValue()):
		ChosenIBIMap.append(4)
		print ChosenIBIMap
	if(self.rb5.GetValue()):
		ChosenIBIMap.append(5)
		print ChosenIBIMap

    def OnAdd(self, event):
        if not self.tc1.GetValue() or not self.tc2.GetValue():
            return
        num_items = self.lc.GetItemCount()
        self.lc.InsertStringItem(num_items, self.tc1.GetValue())
        self.lc.SetStringItem(num_items, 1, self.tc2.GetValue())
	a = int(self.tc1.GetValue())
	b = int(self.tc2.GetValue())
	coords.append([a, b])
	print a*b
        self.tc1.Clear()
        self.tc2.Clear()
	

    def OnRemove(self, event):
        index = self.lc.GetFocusedItem()
        self.lc.DeleteItem(index)


    def OnProcess(self, event):

	
	print coords
	print coords[0]
	print coords[0][1]
	for i in range(len(coords)):
		print 'i= ',i
		indexx=coords[i][0]
		indexy=coords[i][1]
		index=indexx+80*(79-indexy)
		print 'index = ',index
		print 'mapindex = ',mapvec[index]
		print 'mapvec[0] = ',mapvec[0]
		dummyindex=mapvec[index]-1
		print 'dummyindex= ',dummyindex
		pixelvec.append(dummyindex)
	print pixelvec

	print 'Reading in a time series...'
	
	for i in range(len(pixelvec)):
		LongTSPixelIndex=pixelvec[i]

		longvec=[]
		for i in range(2*71990):
			j=2*i*6400+HeaderLength*2+12800*(3)+2*LongTSPixelIndex
			l=2*i*6400+HeaderLength*2+2+12800*(3)+2*LongTSPixelIndex
			datum=struct.unpack("h",self.f1[j:l])[0]
			longvec.append(datum)
		print 'done'
		LongTSPixelMat.append(longvec)
		print 'FilteredLongTS'
		print(FilteredLongTS)
		file = open("nevinTS.txt", "w")
		for i in range(len(longvec)):
			file.write(str(longvec[i])+"\n")
		file.close()
		print 'finished'

	
#######################################
#######FILTER##########################

		dt=0.025
		hpcf = .15
		hpsf = .05
		lpcf = 6
		lpsf = 8


		Nyq = 1/(2*dt)
		Rp = 2
		Rs = 20
		Wp = hpcf/Nyq
		Ws = hpsf/Nyq
		[n,Wn] = buttord(Wp,Ws,Rp,Rs)

		t=range(len(longvec))
		s=longvec-mean(longvec)

		[b,a] = butter(n,Wn,btype='highpass')
		xlp = transpose(lfilter(b,a,transpose(s)))

		Wp = lpcf/Nyq
		Ws = lpsf/Nyq
		[n,Wn] = buttord(Wp,Ws,Rp,Rs)

		[b,a] = butter(n,Wn,btype='lowpass')
		xlpa = transpose(lfilter(b,a,transpose(xlp)))
		FilteredLongTS.append(xlpa)
		

####FILTER######
################
		

    def OnClear(self, event):
        self.lc.DeleteAllItems()
	del coords[:]
	del pixelvec[:]
	del LongTSPixelMat[:]	
	del FilteredLongTS[:]
	#coords =[]
	snapshot = []
	snapvec=[]
	snapvecReorg=[];
	rangevec = []
	snapshotReorg = []
	rangevecReorg = []
	#pixelvec=[]
	#LongTSPixelMat=[]
	#FilteredLongTS=[]



    def OnShowIBIs(self, event):
        frame = CanvasFrameIBIs()
        frame.Show(true)

    def OnShowRetMap(self, event):
        frame = CanvasFrameRetMap()
        frame.Show(true)

    def OnShowTS(self, event):
        frame = CanvasFrameTS()
        frame.Show(true)

    def OnShowTSa(self, event):
        frame = CanvasFrameMap()
        frame.Show(True)

    def OnShowSlider(self, event):
        frame = Slider(None, -1, 'Thresholds')
        frame.Show(True)
        frame.Centre()


    def OnReadData(self,event):
	print self.FileName


	self.f=open(self.FileName,"rb")
	w=open("TS.txt","w");
	self.f1=self.f.read(2000000000)
	print 'done\n'


	f=self.f
	f1=self.f1


	header=struct.unpack(HeaderLength*"h",self.f1[0:HeaderLength*2])
	NumFrames=header[2000] + header[2001]*32000; 
	print NumFrames

#Creating a snapshot: initial vector
	print 'Reading and rearranging an image vector...'
	snapvec=[]
	for i in range(6400):
		j=2*HeaderLength+2*i+3*12800
		l=j+2
		a=struct.unpack("h",self.f1[j:l])[0]
		snapvec.append(a)
	print 'snapvec created'

#Reorganising the snapshot vector to form an image vector
	snapvecReorg=[]
	for i in range(4*QSize**2):
	        a=mapvec[i]-1
	        snapvecReorg.append(snapvec[a])
	print 'done\n'

###ukygdkawufgkaufausdyfgksduygfkuygkasudycfg



	dummy=[[]]
	for i in range(6399):
		dummy.append([])
	print 'dummy recreated done\n'

	print 'Reading 500 frames into the data matrix...'
	for i in range(3200000):
		k=i % 6400
		j=2*i+HeaderLength*2+12800*3
		l=2*i+HeaderLength*2+2+12800*3
		datum=struct.unpack("h",self.f1[j:l])[0]
		dummy[k].append(datum)
	print 'done\n'
#creating a snapshot from the data read in to the data matrix
	print 'Creating snapshot'
	snapshot=[]
	for i in range (6400):
		a=dummy[i][1]
		snapshot.append(a)
	
	print 'Creating vector of pixel ranges'
	rangevec=[]
	for i in range(6400):
		j=max(dummy[i])-min(dummy[i])
		rangevec.append(j)
	
	snapshotReorg=[]
	for i in range(6400):
		a=mapvec[i]-1
		snapshotReorg.append(snapshot[a])
	
	#Rearranging the range vector information to form an image 	vector
#	rangevecReorg=[]
	for i in range(4*QSize**2):
	        a=mapvec[i]-1
	        rangevecReorg.append(rangevec[a])
	
	
	#Normalising the vectors of dummy
	normaldummy=dummy
	for i in range(6400):
		normaldummy[i]=dummy[i]-mean(dummy[i])
	

	#find the max and min values of ranges	
	maxi=max(rangevec)
	for i in range(len(rangevec)):
		if rangevec[i]==maxi:
			maxindex=i
			print 'The max-range is: ',maxi
			print 'An index of the max-range pixel is: ',i
	mini=min(rangevec)
	for i in range(len(rangevec)):
		if rangevec[i]==mini:
			minindex=i
			print 'The min-range is: ',mini
			print 'An index of the min-range pixel is: ',i



    def message(self, event):
        dlg = wx.MessageDialog(self, 'To save one life is as if you have saved the world.', 'Talmud', wx.OK|wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    def choosecolor(self, event):
        dlg = wx.ColourDialog(self)
        dlg.GetColourData().SetChooseFull(True)
        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetColourData()
            self.SetStatusText('You selected: %s\n' % str(data.GetColour().Get()))
        dlg.Destroy()

    def openfile(self, event):
       	dlg = wx.FileDialog(self, "Choose a file", os.getcwd(), "", "*.*", wx.OPEN)
       	if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                mypath = os.path.basename(path)
                self.SetStatusText("Current Data File: %s" % path)
	self.FileName=path
	FileString.append(path)
       	dlg.Destroy()

    def pagesetup(self, event):
        dlg = wx.PageSetupDialog(self)
        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetPageSetupData()
            tl = data.GetMarginTopLeft()
            br = data.GetMarginBottomRight()
            self.SetStatusText('Margins are: %s %s' % (str(tl), str(br)))
        dlg.Destroy()

    def choosefont(self, event):
        default_font = wx.Font(10, wx.SWISS , wx.NORMAL, wx.NORMAL, False, "Verdana")
        data = wx.FontData()
        if sys.platform == 'win32':
            data.EnableEffects(True)
        data.SetAllowSymbols(False)
        data.SetInitialFont(default_font)
        data.SetRange(10, 30)
        dlg = wx.FontDialog(self, data)
        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetFontData()
            font = data.GetChosenFont()
            color = data.GetColour()
            text = 'Face: %s, Size: %d, Color: %s' % (font.GetFaceName(), font.GetPointSize(),  color.Get())
            self.SetStatusText(text)
        dlg.Destroy()

    def opendir(self, event):
        dlg = wx.DirDialog(self, "Choose a directory:", style=wx.DD_DEFAULT_STYLE | wx.DD_NEW_DIR_BUTTON)
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You selected: %s\n' % dlg.GetPath())
        dlg.Destroy()

    def singlechoice(self, event):
        sins = ['Greed', 'Lust', 'Gluttony', 'Pride', 'Sloth', 'Envy', 'Wrath']
        dlg = wx.SingleChoiceDialog(self, 'Seven deadly sins', 'Which one?', sins, wx.CHOICEDLG_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You chose: %s\n' % dlg.GetStringSelection())
        dlg.Destroy()

    def textentry(self, event):
        dlg = wx.TextEntryDialog(self, 'Enter some text','Text Entry')
        dlg.SetValue("Default")
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You entered: %s\n' % dlg.GetValue())
        dlg.Destroy()

class MyApp(wx.App):
    def OnInit(self):
        frame = MyFrame(None, -1, 'Analysis of RS Data')
        frame.Show(True)
        frame.Centre()
        return True

app = MyApp(0)
app.MainLoop()
