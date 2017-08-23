from __future__ import division

import numpy as np
import ctypes
import os

import sys
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLineEdit, QGridLayout, QToolTip, QPushButton, QSlider
from PyQt5.QtGui import QIcon
import PyQt5
from matplotlib.backends.backend_qt5agg import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import time
import threading

''' 
Author - Daniel J. Whiting 
Date modified - 10/08/2017
--- Installation ---
Requires standard 64-bit python 2 distribution (e.g. anaconda) including PyQt5 library
Install 64 bit thorcam package
Change PATH variable to point to installation location of thorcam
--- Usage ---
Currently continuously reads and displays images from camera with adjustable exposure
Clicking button to calculate waist will print result to terminal window
--- Developer notes ---
Functions return 0 if successful or -1 if failed, 125 means invalid parameter
--- Changelog ---
Defined new function 'Gaussian' to take place of lambda functions
Corrected gaussian formula r->2r
Added horizontal and vertical slices, centered on the peak intensity
Added display of last 20 waists recorded.
Added toggle to show/hide the above features. 
'''




class IS_RECT(ctypes.Structure):
    _fields_ = [
        ("s32X", ctypes.c_int),
        ("s32Y", ctypes.c_int),
        ("s32Width", ctypes.c_int),
        ("s32Height", ctypes.c_int)
    ]


class cameraAPI():
    def __init__(self):
        # Load DLL into memory
        PATH = r'C:\Program Files\Thorlabs\Scientific Imaging\ThorCam'
        os.environ['PATH'] = ';'.join([PATH, os.environ['PATH']])
        self.dll = ctypes.CDLL(os.path.join(PATH, 'uc480_64.dll'))

        # number_of_cameras = ctypes.c_int(0)
        # _chk(self.clib.is_GetNumberOfCameras(byref(number_of_cameras)))
        # if number_of_cameras.value < 1:
        #    raise RuntimeError("No camera detected!")

        self.ModuleHandle = ctypes.c_int()
        self.dll.is_InitCamera(ctypes.pointer(self.ModuleHandle))

        # Set AOI
        rectAOI = IS_RECT()
        self.dll.is_AOI(self.ModuleHandle, 2, ctypes.pointer(rectAOI), 4 * 4)
        self.shape = (rectAOI.s32Width, rectAOI.s32Height)

        # Setting monocrome 8 bit color mode
        self.dll.is_SetColorMode(self.ModuleHandle, 6)

        # Allocate memory:
        self.pid = ctypes.c_int()
        self.ppcImgMem = ctypes.c_char_p()
        self.dll.is_AllocImageMem(self.ModuleHandle, self.shape[0], self.shape[1], 8, ctypes.pointer(self.ppcImgMem),
                                  ctypes.pointer(self.pid))  # 10 bit not 8?
        self.dll.is_SetImageMem(self.ModuleHandle, self.ppcImgMem, self.pid)

        # Set exposure
        # self.get_exposure_range()
        # self.update_exposure_time(0.5) # in ms

        # Other settings
        self.dll.is_SetExternalTrigger(self.ModuleHandle, 8)
        self.dll.is_SetHardwareGain(self.ModuleHandle, 0, 0, 0, 0)
        self.dll.is_EnableAutoExit(self.ModuleHandle, 1)

    def get_exposure_range(self):
        self.exposure_min = ctypes.c_double()
        self.exposure_max = ctypes.c_double()
        self.exposure_step = ctypes.c_double()
        self.dll.is_Exposure(self.ModuleHandle, 8, ctypes.pointer(self.exposure_min), 8)
        self.dll.is_Exposure(self.ModuleHandle, 9, ctypes.pointer(self.exposure_max), 8)
        self.dll.is_Exposure(self.ModuleHandle, 10, ctypes.pointer(self.exposure_step), 8)
        self.exposure_min = self.exposure_min.value
        self.exposure_max = self.exposure_max.value
        self.exposure_step = self.exposure_step.value
        self.exposure_min = 0.037
        self.exposure_max = 983
        self.exposure_step = 0.037

    def update_exposure_time(self, t, units='ms'):
        """Set the exposure time."""
        IS_EXPOSURE_CMD_SET_EXPOSURE = 12
        nCommand = IS_EXPOSURE_CMD_SET_EXPOSURE
        Param = ctypes.c_double(t)
        SizeOfParam = 8
        self.dll.is_Exposure(self.ModuleHandle, nCommand, ctypes.pointer(Param), SizeOfParam)

    # def toggle_auto_exposure(self):
    #			#par1 = ctypes.c_double(0) # ctypes.pointer(par1)
    # self.dll.is_SetAutoParameter(self.ModuleHandle,0x8802,ctypes.pointer(par1),0) # Switch off auto expose

    def get_image(self):
        # Allocate memory for image:
        img_size = self.shape[0] * self.shape[1]
        c_array = ctypes.c_char * img_size
        c_img = c_array()

        # Take one picture: wait time is waittime * 10 ms:
        waittime = ctypes.c_int(100)
        self.dll.is_FreezeVideo(self.ModuleHandle, waittime)

        # Copy image data from the driver allocated memory to the memory that we allocated.
        self.dll.is_CopyImageMem(self.ModuleHandle, self.ppcImgMem, self.pid, c_img)

        # Convert to python array
        img_array = np.frombuffer(c_img, dtype=ctypes.c_ubyte)
        img_array.shape = (self.shape[1], self.shape[0])
        return img_array

        # def save_image(self):


class main(QWidget):
	def __init__(self):
		QWidget.__init__(self)

		self.cam = cameraAPI()

		self.activateExtra = 1
		self.continuous = 0

		self.waistListX = np.zeros(20)
		self.waistListY = np.zeros(20)

		# self.setGeometry(300, 300, 300, 220)
		self.setWindowTitle('Camera Software')
		self.setWindowIcon(QIcon('web.png'))

		# Add Matplotlib Canvas to plot ThorCam image to.
		self.fig = plt.figure(figsize=(5, 5))
		self.canvas = FigureCanvasQTAgg(self.fig)
		# Add second canvas to print the histogram to.
		self.intensityFig = plt.figure(figsize=(5, 2))
		self.intensityCanvas = FigureCanvasQTAgg(self.intensityFig)

		self.waistTextBox = QLineEdit(self)

		Button_1 = QPushButton('Calc waist', self)
		Button_1.clicked.connect(self.calc_waists)

		self.buttonContinous = QPushButton('Toggle Continuous Mode', self)
		self.buttonContinous.clicked.connect(self.toggleContinuousMode)

		ButtonShowHide = QPushButton('Toggle Graphs', self)
		ButtonShowHide.clicked.connect(self.showHide)

		self.Exposure_slider = QSlider(orientation=Qt.Horizontal, parent=self)
		self.Exposure_slider.setMinimum(1)
		self.Exposure_slider.setMaximum(100)
		self.Exposure_slider.setValue(50)
		self.On_exposure_change()
		self.Exposure_slider.valueChanged.connect(self.On_exposure_change)

		# set the layout
		layout = QGridLayout()
		# layout.addWidget(self.toolbar)
		layout.addWidget(Button_1, 1, 0)
		layout.addWidget(self.Exposure_slider, 2, 0)
		layout.addWidget(self.canvas, 3, 0)
		layout.addWidget(ButtonShowHide, 4, 0)
		layout.addWidget(self.waistTextBox, 2, 1)
		layout.addWidget(self.intensityCanvas, 3, 1)
		layout.addWidget(self.buttonContinous,4,1)
		# layout.addWidget(self.button)

		layout.addWidget

		self.showHide()

		self.setLayout(layout)

		self.show()

		self.run_stream = True
		self.camera_stream()

	def toggleContinuousMode(self):
		self.continuous = (self.continuous + 1) % 2
		return

	def showHide(self):
		if self.activateExtra == 1:
			self.intensityCanvas.hide()
			self.waistTextBox.hide()
			self.buttonContinous.hide()
		else:
			self.intensityCanvas.show()
			self.waistTextBox.show()
			self.buttonContinous.show()
		self.activateExtra = (self.activateExtra + 1) % 2

	def On_exposure_change(self):
		new_exposure = 0.037 * 10 ** (self.Exposure_slider.value() / 23)
		self.cam.update_exposure_time(new_exposure)

	def closeEvent(self, event):
		self.run_stream = False
		time.sleep(1)
		event.accept()  # let the window close

	def camera_stream(self):
		t = threading.Timer(0, function=self.capture_image)
		t.daemon = True
		t.start()

	def get1DIntensity(self, axis):

		maxIndex = np.argmax(self.imdata)

		maxXIndex, maxYIndex = np.unravel_index(maxIndex, self.imdata.shape)

		if axis == 'v':
			oneDIntensity = self.imdata[maxXIndex, :]
		if axis == 'h':
			oneDIntensity = self.imdata[:, maxYIndex]

		return oneDIntensity, (maxXIndex, maxYIndex)

	def capture_image(self):
		# Create the matplotlib axis to display the image data
		self.ax = self.fig.add_subplot(111)

		from scipy import misc
		self.imdata = self.cam.get_image()
		#from scipy.misc import imread
		#self.imdata = imread('./gaussian.jpg')

		self.image = self.ax.imshow(self.imdata, vmin=0, vmax=255, cmap='gray', origin='lower left')

		# Create the matplotlib axis to display the histogram of intensities
		self.hax = self.intensityFig.add_subplot(311)
		self.hax.set_title('Horizontal Histogram')
		self.hdata = self.get1DIntensity('h')[0]
		self.hplot, = self.hax.plot(self.hdata, color = '0.8')
		self.hintplot, = self.hax.plot(np.zeros(np.sum(self.imdata, axis=1).shape), color = '0.5', linewidth = 3)
		self.hintfit, = self.hax.plot(np.zeros(np.sum(self.imdata, axis=1).shape), color = 'g')
		self.hax.set_ylim(0,255)

		
		self.vax = self.intensityFig.add_subplot(312)
		self.vax.set_title('Vertical Histogram')
		self.vdata = self.get1DIntensity('v')[0]
		self.vplot, = self.vax.plot(self.vdata, color = '0.8')
		self.vintplot, = self.vax.plot(np.sum(self.imdata, axis=0), color = '0.5', linewidth = 3)
		self.vintfit, = self.vax.plot(np.zeros(np.sum(self.imdata, axis=0).shape), color = 'g')
		self.vax.set_ylim(0,255)

		# Create axis to display waists
		self.wax = self.intensityFig.add_subplot(313)
		self.wax.set_title('Last 20 Waists')

		self.wxplot, = self.wax.plot(self.waistListX)
		self.wyplot, = self.wax.plot(self.waistListY)

		while self.run_stream:
			self.imdata = self.cam.get_image()
			self.hdata, self.hmax = self.get1DIntensity('h')
			self.vdata, self.vmax = self.get1DIntensity('v')
			self.image.set_data(self.imdata)
			
			self.hplot.set_ydata(self.hdata)
			self.vplot.set_ydata(self.vdata)
			
			vint = np.sum(self.imdata, axis=0).astype(np.float64)
			hint = np.sum(self.imdata, axis=1).astype(np.float64)
			
			vint *= 255/vint.max()
			hint *= 255/hint.max()
			
			self.hintplot.set_ydata(hint)
			self.vintplot.set_ydata(vint)
			
			self.wxplot.set_ydata(self.waistListX)
			self.wyplot.set_ydata(self.waistListY)


			if self.continuous ==1:
				self.calc_waists()

			self.wax.relim()
			self.wax.autoscale_view(True, True, True)


			self.canvas.draw()
			self.intensityCanvas.draw()

			self.canvas.flush_events()
			self.intensityCanvas.flush_events()
			self.intensityFig.tight_layout()

	def gaussian(self,x, a, x0, b, wx):
		a = np.abs(a)
		return a * np.exp(-2*((x - x0) / wx) ** 2)+ b

	def calc_waists(self):
		try:
			xdata = np.sum(self.imdata, axis=1)  # x
			ydata = np.sum(self.imdata, axis=0)  # x
			
			xaxis = np.arange(len(xdata))
			yaxis = np.arange(len(ydata))

			p0x = (xdata.max(), 600, 0, 600)
			p0y = (ydata.max(), 600, 0, 600)

			px, covx = curve_fit(self.gaussian, xaxis, xdata,
						   p0=p0x)
			py, covy = curve_fit(self.gaussian, yaxis, ydata,
						   p0=p0y)
						 
			hfit = self.gaussian(xaxis, *px)
			vfit = self.gaussian(yaxis, *py)
						 
			hfit *= 255./hfit.max()
			vfit *= 255./vfit.max()
			
			self.hintfit.set_ydata(hfit)
			self.vintfit.set_ydata(vfit)

			wx = np.abs(px[-1])
			wy = np.abs(py[-1])


			pixel_size = 5.2e-3  # mm

			self.waistListX = np.roll(self.waistListX, 1)
			self.waistListY = np.roll(self.waistListY, 1)
			self.waistListX[0] = wx * pixel_size
			self.waistListY[0] = wy * pixel_size

			message = 'wx = ' + str(wx * pixel_size) + ' | wy = ' + str(wy * pixel_size) + ' (mm)'
			self.waistTextBox.setText(message)
			print(message)
		except Exception as e:
			print( e )
			None


if __name__ == '__main__':
	app = QApplication(sys.argv)
	main = main()
sys.exit(app.exec_())