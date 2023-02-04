% Copyright (C) 2017 - .... Burkhard Schmidt 
%               2009 - 2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init(delay, channel)
global atomic hamilt space time

prt.disp ('****************************************')
prt.disp ('Calculation of the pump probe signal for')
prt.disp ('excitation of Na2.                      ')
prt.disp ('see J.Chem.Phys. 100:5448               ')
prt.disp ('****************************************')

% Note that the function argument "delay" gives the time delay between
% pump and probe pulse in units of 10 fs. This affects the calculation 
% in a number of ways:
%
% 1. We place the 30 fs probe pulse at the given time delay.
% 2. We propagate from 100 fs before the probe to 100 fs after probe.
% 3. We take the initial wavefunction from the pump-only calculation.
%
% 2. and 3. speed up the calculation a _lot_. (Delays are up to 2 ps!)
probe_time = delay + 10;    % Necessary for some index computations. Takes care
                            % of the fact that the pump pulse is centered at
                            % 100 fs.

% Number of coupled equations
hamilt.coupling.n_eqs     = 3;
hamilt.coupling.represent = 'dia';
hamilt.coupling.labels    = {'X^1\Sigma_g^+', 'A^1\Sigma_u^+', '2^1\Pi_g'};
hamilt.coupling.ini_rep   = 'dia';
hamilt.coupling.ini_norm  = false;       % don't normalize initial wave function
                                         % not sure why this is needed here
% Grid definition
space.dof{1}       = dof.fft;            % using FFT grid
space.dof{1}.n_pts = 256;                % number of points
space.dof{1}.x_min = 4;                  % lower boundary of the grid
space.dof{1}.x_max = 14;                 % upper boundary of the grid
space.dof{1}.mass  = atomic.mass.Na23/2; % reduced mass of 23Na2

% Temporal discretisation
time.steps.m_start  = probe_time - 10;   % index of first time step
time.steps.m_stop   = probe_time + 10;   % index of last time step
time.steps.m_delta  = 10/atomic.t.fs;    % 10 fs per time step
time.steps.s_number = 500;               % propagation time step 20 as

% Two identically shaped Gaussian pulses
time.pulse{1}       = efi.gauss;         % Gaussian-shaped pulse
time.pulse{1}.fwhm  = 30/atomic.t.fs;    % with 30 fs FWHM
time.pulse{1}.delay = 100/atomic.t.fs;   % delayed by 100 fs
time.pulse{1}.phase = 0;                 % locked phase
time.pulse{1}.frequ = 0.073;             % frequency (corr. to 625 nm)
time.pulse{1}.ampli = sqrt(10/atomic.I.GW_cm2); % 10 GW/cm^2

time.pulse{2}       = efi.gauss;         % Gaussian-shaped pulse
time.pulse{2}.fwhm  = 30/atomic.t.fs;    % with 30 fs FWHM
time.pulse{2}.delay = 10*probe_time/atomic.t.fs;  % delayed by 100 fs + "delay"
time.pulse{2}.phase = 0;                 % locked phase
time.pulse{2}.frequ = 0.073;             % frequency (corr. to 625 nm)
time.pulse{2}.ampli = sqrt(10/atomic.I.GW_cm2); % 10 GW/cm^2 

% Hamiltonian operator
hamilt.truncate.e_min = -0.03;           % lower truncation of energy
hamilt.truncate.e_max =  0.25;           % upper truncation of energy

for m=1:3
    hamilt.pot{m,m}          = pot.interp;   % interpolate tabulated potential
    hamilt.pot{m,m}.pos_conv = atomic.l.A;   % conversion factor for coordinates
    hamilt.pot{m,m}.pot_conv = atomic.E.eV;  % conversion factor for energies
end

hamilt.dip{1}             = cell(3);
hamilt.dip{1}{1,2}        = dip.taylor;  % transition dipole moment
hamilt.dip{1}{1,2}.vshift = 1;           % constant ONE
hamilt.dip{1}{2,3}        = dip.taylor;  % transition dipole moment
hamilt.dip{1}{2,3}.vshift = 1;           % constant ONE

% Initial wave function; Load from pump-only run.
time.corr         = init.load;           % load from file
time.corr.dir     = '../2/';
time.corr.file    = 'pump';
time.corr.index   = delay + 1;      % index 1 is the initial wave function!
time.corr.channel = channel;

