% Copyright (C) 2008-2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

%% obtain the spectrum from the autocorrelation of a TDSE solution

% Run the calculation
qm_setup();
qm_init();
qm_propa();

% Second part: Do something clever with the autocorrelation function.

% First: Window the autocorrelation function, which will smooth out the
% peaks in the Fourier transform to Gaussians of width ~ 1/tau
global atomic time

tau = 100000;

autocorr = time.steps.acf .* exp(-(time.steps.s_grid/tau).^2);

% I'll never get the hang of discrete Fourier transforms. Anyway, all we
% need to do is to cut of the last element of the grids, then we have even
% numbers and can just copy the code from the fft grid initialisation.

autocorr = autocorr(1:end-1);
ftrafo   = real(ifft(autocorr, [], 1));

minfreq = -pi / (time.steps.s_grid(2) - time.steps.s_grid(1));
maxfreq = -minfreq;
frequencies = linspace(minfreq, maxfreq - 2*maxfreq/(length(autocorr)), length(autocorr));
frequencies = fftshift(frequencies);

% ftrafo contains the smoothed Fourier-transform of the autocorrelation, and
% frequencies the corresponding frequency values (== energies) in atomic units.
% Due to the internals of fft, they are output from w=0 to w_max, and w_min to w_-1.

% For easier comparison, convert the energy scale to cm^-1
global atomic
frequencies = frequencies * atomic.w.cm_1;

% Now we are only interested in properties between 5000 and 25000 cm^-1, so remove
% everything outside.
ftrafo(frequencies > 25000) = 0;
ftrafo(frequencies < 5000) = 0;
frequencies(frequencies > 25000) = 25001;
frequencies(frequencies < 5000) = 4999;

% Plot the result
clf;
plot(frequencies, ftrafo, 'LineWidth', 2);
saveas(gcf, 'autocorr.jpg');

qm_cleanup();

%% obtain the spectrum from direct solution of the TISE
qm_setup();
qm_init();
qm_bound();
qm_cleanup();

global expect
format long
sort(expect.total*atomic.w.cm_1)
