import requests
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
import os
import h5py
from pi_stage import PI_stage

from small_lab_gui.helper import bokeh_gui_helper as bgh
from small_lab_gui.helper import bokeh_plot_helper as bph

longstage = False

if longstage:
    from smc100pp import SMC100
else:
    from pi_stage import PI_stage

class osc_alignment():
    def __init__(self):
        self.logbook = logbook
        self.title = 'Alignment'
        self.totalCount = 0


class oscilloScan():
    def __init__(self, delayer, osc) -> None:

        self.logbook = logbook
        self.title = 'Measurement'

        self.delayer = delayer
        self.osc = osc

        # pcolor plot to display results
        self.imdata = 0
        self.sumImageSpectrum = 0
        self.sumSpectrum = 0
        self.totalFFT = 0
        self.lastSpectrum = 0
        self.lastFFT = 0
        self.allSpectrum = 0
        self.scanNum = 0

        self.imagePlot = bph.plot_false_color()
        self.imagePlot.image()
        #sum of each scan
        self.sumImagePlot = bph.plot_2d()
        self.sumImagePlot.line()
        #add togethor
        self.sumPlot = bph.plot_2d()
        self.sumPlot.line()
        #FFT of total scan
        self.sumFFTPlot = bph.plot_2d()
        self.sumFFTPlot.line()
        #last scan
        self.lastPlot = bph.plot_2d()
        self.lastPlot.line()
        #FFT of last scan
        self.lastFFTPlot = bph.plot_2d()
        self.lastFFTPlot.line()

        #buttons and layout
        self.startBtn = bokeh.models.widgets.Button(
            label='Start', button_type='success')
        self.stageStartPos = bokeh.models.widgets.TextInput(
                title='Start Position [um]', value='50')        
        #experiment parameter
        self.sweepTimeInput = bokeh.models.widgets.TextInput(
            title='Integration time [s]', value='100') #seconds
        self.delay = np.arange(10000)*(100/10000)*2*3.33564*10**-15
        self.numCycle = 1 # number of cycle after initiation of the wave generator
        self.numSweeps = bokeh.models.widgets.TextInput(
            title='Sweep Number', value='30') # number of sweeps

        self.inputs = bokeh.layouts.Column(
            self.startBtn, self.sweepTimeInput, self.numSweeps, self.stageStartPos
        ) 
        
        self.layout = bokeh.layouts.row(
            self.inputs,
            bokeh.layouts.column(
                bokeh.layouts.row(self.lastPlot.element, self.lastFFTPlot.element, self.imagePlot.element),
                bokeh.layouts.row(self.sumPlot.element, self.sumFFTPlot.element, self.sumImagePlot.element)
            ),
            width=800)

    def start(self):
        #open shutter
        #self.shutter.set_shutter_mode(modes=['F','O'])

        # in case this measurement is running, stop it
        if self.running.am_i_running(self):
            self.stop()
        # in case a different measurement is running, do nothing
        elif self.running.is_running():
            pass
        else:
            self.running.now_running(self)
            self.startBtn.label = 'Stop'
            self.startBtn.button_type = 'danger'

            #in future here should be able to set the phase of the wave output of wave generator
            self.delayer.setSweepTime(int(self.sweepTimeInput.value))
            self.delayer.setCycleNum(int(self.numCycle))#number of the cycle after initiate the wave generator
            
            #oscilloscope
            self.osc.checkDataSource()
            #in the future here should be able to remotely set the parameter of the osciloscope by the IP address.

            # scan start time for save name
            self.now = datetime.datetime.now()

            self.imdata = 0
            self.sumSpectrum = 0
            self.totalFFT = 0
            self.lastSpectrum = 0
            self.lastFFT = 0
            self.allSpectrum = 0
            self.scanNum = 0

            # switch start to stop button
            self.startBtn.label = 'Stop'
            self.startBtn.button_type = 'danger'

            # create the measurment thread
            self.measurement = measurement.measurement(
                inputs=None,
                sequence=[
                    self.delayer.initWaveGen,
                    measurement.sleep_function(int(self.sweepTimeInput.value)+2),
                    self.osc.requestData,
                    self.delayer.move_absolute_um(float(self.startPos.value)),
                    measurement.sleep_function(0.3)
                    ],
                update=measurement.bokeh_update_function(
                    self.update, self.doc),
                init=self.delayer.move_absolute_um(float(self.startPos.value)),
                finish=self.delayer.close(),
                save_output=False)

    def update(self, data):
        self.im_data = np.array([d[1][0] for d in data])
        self.lastSpectrum = self.im_data[-1]
        plot_data = np.transpose(im_data)
        self.sumSpectrum = np.sum(plot_data)
        self.sumImageSpectrum = np.sum(plot_data, 1)

        self.scanNum += 1

        # update plots
        try:
            self.lastPlot.update(
                num=0, x=self.delay, y=self.lastSpectrum
            )
            self.imagePlot.update(
                num=0, x=self.delay*1e15, y=np.arange(self.scanNum), z=plot_data
            )
            self.sumImagePlot.update(
                num=0, x=self.delay*1e15, y=self.sumSpectrum
            )
        except:
            print("plot error!")

class oscScope():
    def __init__(self, dataSource = 'CH1', IP = '129.27.156.205') -> None:
        self.dataSource = dataSource
        self.IP = IP
        self.commandIP = 'http://129.27.156.205/Comm.html'
    
    def checkDataSource(self):
        if dataSource == self.dataSource:
            print('The data source is <' + self.dataSource + '>.')
            pass
        else:
            response = requests.post(
                'http://129.27.156.205/Comm.html', data={'COMMAND': 'DATA:SOURCE '+self.dataSource})
            print('The data source is modified to <' + self.dataSource + '>.')
    
    def requestData(self):
        response = requests.post(
            'http://129.27.156.205/Comm.html', data={'COMMAND': 'CURVe?'})
        r = response.text
        endOfHeader = r.find("NAME=\"name\""+">")+len("NAME=\"name\""+">")
        endOfData = r.find("</TEXTAREA>")
        data = r[endOfHeader:endOfData]  # remove the header info from the response
        data = np.fromstring(data, sep=',')  # 10k points
        return data

class shutter_gui():
    def __init__(self, doc):
        self.shutter = Shutter('COM25')
        self.title = 'Shutter'
        # bokeh doc for callback
        self.doc = doc

        self.width = 300
        self.height = round(0.618*self.width)
        self.openBtn = bokeh.models.widgets.Button(
            label='Open All', button_type='success', width=self.width, height=self.height)
        self.openFiberBtn = bokeh.models.widgets.Button(
            label='Open Fiber', button_type='success', width=self.width, height=self.height)
        # start thread callback
        self.openBtn.on_click(self.openAll)
        self.openFiberBtn.on_click(self.openFiber)

        self.layout = bokeh.layouts.row(self.openBtn, self.openFiberBtn, width=1600, height=400)

    
    def openAll(self):
        if self.openBtn.label == 'Open All':
            self.shutter.set_shutter_mode(modes=['F', 'O'])
            # switch start to stop button
            self.openBtn.label = 'Close All'
            self.openBtn.button_type = 'danger'
            self.openFiberBtn.button_type = 'danger'
        else:
            self.closeAll()

    def openFiber(self):
        if self.openFiberBtn.label == 'Open Fiber':
            self.shutter.set_shutter_mode(mode='F')
            self.openFiberBtn.label = 'Close Fiber'
            self.openFiberBtn.button_type = 'danger'
            self.openBtn.button_type = 'success'
        else:
            self.closeFiber()

    def closeAll(self):
        self.shutter.set_shutter_mode(modes=['K', 'P'])
        self.openBtn.label = 'Open All'
        self.openBtn.button_type = 'success'
        self.openFiberBtn.button_type = 'success'


    def closeFiber(self):
        self.shutter.set_shutter_mode(mode='K')
        self.openFiberBtn.label = 'Open Fiber'
        self.openFiberBtn.button_type = 'success'
        self.openBtn.button_type = 'success'

    def close(self):
        self.closeAll()

class osc_session_handler(bgh.bokeh_gui_session_handler):
    def open_session(self, doc):
        self.running = bgh.running()

        # hardware
        osc = oscScope()

        if longstage:
            stage = SMC100(1, 'COM22', silent=True)
        else:
            stage = PI_stage('COM23')

        # open logbook to auto-save scans
        # logbook = elog(host='localhost', port=8080, logbook='demo')



        shuttergui = shutter_gui(doc=doc)

        self.title = 'TOF Readout'
        self.tabs = [
            {'layout': alignmentgui.layout, 'title': alignmentgui.title},
            {'layout': measurementgui.layout, 'title': measurementgui.title},
            {'layout': shuttergui.layout, 'title': shuttergui.title}]

        # this list is auto-closed, all close functions of the
        # added objects are called at session destruction
        self.close_list.append(measurementgui)
        self.close_list.append(shuttergui)
        self.close_list.append(stage)



numCycle = 1 # number of cycle after initiation of the wave generator
numSweeps = 40 # number of sweeps
stage = PI_stage('COM23')


now = datetime.datetime.now()
os.makedirs(
    now.strftime('%Y-%m') + '/'
    + now.strftime('%Y-%m-%d'), exist_ok=True)
fname = (now.strftime('%Y-%m') + '/'
         + now.strftime('%Y-%m-%d')
         + '/scan_osc_'
         + now.strftime('%Y-%m-%d-%H-%M-%S'))
with h5py.File(fname + '.hdf5', 'w') as f:
    for i in range(numSweeps):
        stage.move_absolute_um(0)
        time.sleep(0.5)
        stage.initWaveGen()
        print('This is the No.' + str(i) + ' sweep.')
        time.sleep(sweepTime+1)
        stage.stopWaveGen()
        response = requests.post(
            'http://129.27.156.205/Comm.html', data={'COMMAND': 'CURVe?'})
        r = response.text
        endOfHeader = r.find("NAME=\"name\""+">")+len("NAME=\"name\""+">")
        endOfData = r.find("</TEXTAREA>")
        data = r[endOfHeader:endOfData]  # remove the header info from the response
        data = np.fromstring(data, sep=',')  # 10k points
        #print(data.shape)

        #plt.plot(data)
        #plt.show(block=False)

        try:
            # save data hdf5

            f.create_dataset('scan'+str(i), data=data)
            f.flush()
        except Exception as e:
            print('save error')
            print(e)

