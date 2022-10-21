%function fre=demodulator(signal, fc, fs, fdev)
%signal = h5read('C:\Users\user\Desktop\Data_newTOF\2021-10\2021-10-10\scan_tof_2021-10-10-22-28-56.hdf5','/daq_readout_1');
signal = dlmread('C:\Users\user\Desktop\Data_newTOF\2022-03\2022-03-27\delayMonitor.dat');
signal = signal(1:300000);

fc=63;
fdev=5;
fs=300;
filter_signal = signal;
fre = pmdemod(filter_signal,fc,fs,fdev);
t=[0:1:(size(fre)-1)]./fs;
plot(t, fre)
%plot(t,[filter_signal*30,bandpass(fre,[58 68],300)])
