% Copyright (C) 2009-2010 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global atomic hamilt plots space time

prt.disp('*********************************************')
prt.disp('Control of vibrational modes of HF via       ')
prt.disp('chirped infrared laser pulses.               ')
prt.disp('                                             ')
prt.disp('Figure 1 in the reference                    ')
prt.disp('                                             ')
prt.disp('see S. Chelkowski, A.D. Bandrauk, P.B. Corkum')
prt.disp('Phys. Rev. Lett. 65, page 2355               ')
prt.disp('*********************************************')

% Reduced mass of 1H-19F
m_F = atomic.mass.F19;
m_H = atomic.mass.H1;
m_r = m_F * m_H / (m_F+m_H);

% Define the grid
space.dof{1}       = dof.fft;            % use an FFT grid
space.dof{1}.mass  = m_r;                % reduced mass 
space.dof{1}.n_pts = 2048;               % number of points
space.dof{1}.x_min = 0.5;                % lower boundary; not given in paper
space.dof{1}.x_max = 130.5;              % upper boundary

% Simple Morse potential
hamilt.pot{1,1}     = pot.morse;
hamilt.pot{1,1}.d_e = 6.125 / atomic.E.eV; % dissociation energy
hamilt.pot{1,1}.r_e = 1.7329;            % equilibrium distance
hamilt.pot{1,1}.alf = 1.1741;            % range parameter

% Linear dipole moment
hamilt.dip{1}{1,1} = dip.interp;         % load from an external file

% Hamiltonian: truncation of energies
hamilt.truncate.e_min = 0;
hamilt.truncate.e_max = 1;

% time step parameters
time.steps.m_start  = 0;                 % index of first time step
time.steps.m_stop   = 120;               % index of last time step
time.steps.m_delta  = 8.41/atomic.t.fs;  % laser cycles (8.41 fs)
time.steps.s_number = 1000;              % number of sub-steps per main step

% Initial wave function is the ground state of the Morse oscillator itself
time.dof{1}     = init.morse;
time.dof{1}.m_r = space.dof{1}.mass;     % reduced mass
time.dof{1}.r_e = hamilt.pot{1,1}.r_e;   % equilibrium distance
time.dof{1}.d_e = hamilt.pot{1,1}.d_e;   % dissociation energy
time.dof{1}.alf = hamilt.pot{1,1}.alf;   % range parameter
time.dof{1}.n_q = 0;                     % we want the groundstate

% Plots: Draw only the wave function, not the Wigner-transform (takes ages)
plots.density       = vis.curve;
plots.density.range = true;              % custom plotting ranges
plots.density.x_min = 0.5;
plots.density.x_max = 6;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 0.2;

% the electric field is composed of four pulses. The turn-on pulse is
% interpolated from a file. After that, we have a plateau, and chirp the
% frequency up to time t_N. We continue for some time without chirping the
% pulse, and the final turn-off pulse is again interpolated from file.
% A major complication in the chirped part of the pulse is the fact that the
% reference expresses the chirp as a function of (t - t_1), while the WavePacket
% implementation wants to have the chirp as a function of (t - delay time), and
% also multiplies the time-dependent frequency with this shifted time,
% where the delay time is the center between t_0 (when the chirping starts) and
% t_8 (when chirping ends). Also, the phases of the electric fields are always
% defined relative to the time delay, so static phase shifts have to be added to
% avoid jumps in the electric field when switching from one field to the next.

% some values we need for the calculations; They can be extracted from the
% reference.
% A  = -0.5224;
B  = 0.0419;
Q  = 7316;
t0 = 9513;
t1 = 12286;
t8 = t1 + 2 * Q * (sqrt(9) - sqrt(2));   % eq. 12 in the reference
w0 = 2 * B * hamilt.pot{1,1}.d_e;        % harmonic frequency of the Morse potential
wi = w0 * (1 - B);                       % initial frequency of the pulse
wf = w0 - w0 * B * ( (t8-t1)^2 / (4*Q^2) + sqrt(2) * (t8-t1) / Q + 3/2);
                                         % final frequency of the pulse

T = (t0 + t8) /2;                        % delay time of the chirped pulse
quadratic = -B*w0 / (2*Q^2);             % quadratic chirp
linear    = -B*w0 / (4*Q^2) * (3*T - 2*t1) - sqrt(2)*B*w0 / Q;
                                         % linear chirp
constant  = -B*w0 / (4*Q^2) * (T^2 - t1^2 + 2*(T-t1)^2) ...
            - sqrt(2)*B*w0 / Q * (2*T - t1) + (w0 - 3/2*B*w0);
                                         % time-independent part of chirped frequency
% static phase of the second and third pulse to ensure a continuous field.
phi2 = -B * w0 / (4*Q^2) * T * (T-t1)^2 ...
        - sqrt(2) * B * w0 / Q * T * (T - t1) + (w0 - B * w0*3/2) * T;
phi3 = wf* t8 + wf*t0/10;

time.pulse{1}           = efi.inter;     % pulse shapes
time.pulse{2}           = efi.recta;     % pulse shapes
time.pulse{3}           = efi.recta;     % pulse shapes
time.pulse{4}           = efi.inter;     % pulse shapes

time.pulse{1}.file      = 'turn_on.dat'; % files to load from
time.pulse{4}.file      = 'turn_off.dat';% files to load from
time.pulse{1}.method    = 'spline';      % interpolation method
time.pulse{4}.method    = 'spline';      % interpolation method

time.pulse{1}.ampli     = sqrt(10/atomic.I.TW_cm2);
time.pulse{2}.ampli     = sqrt(10/atomic.I.TW_cm2);
time.pulse{3}.ampli     = sqrt(10/atomic.I.TW_cm2);
time.pulse{4}.ampli     = sqrt(10/atomic.I.TW_cm2);

time.pulse{1}.delay     = t0;            % delay time
time.pulse{2}.delay     = T;             % delay time
time.pulse{3}.delay     = t8+t0/10;      % delay time
time.pulse{4}.delay     = t8+t0/5;       % delay time

time.pulse{1}.frequ     = wi;            % central frequency
time.pulse{2}.frequ     = constant;      % central frequency
time.pulse{3}.frequ     = wf;            % central frequency
time.pulse{4}.frequ     = wf;            % central frequency

time.pulse{2}.linear    = linear;        % linear chirp
time.pulse{2}.quadratic = quadratic;     % quadratic chirp

time.pulse{1}.phase     = wi*t0;         % match phases
time.pulse{2}.phase     = phi2;          % match phases
time.pulse{3}.phase     = phi3;          % match phases
time.pulse{4}.phase     = phi3+wf*t0/10; % match phases

time.pulse{1}.fwhm      = 0;             % FWHM (/length) of pulses
time.pulse{2}.fwhm      = t8-t0;         % FWHM (/length) of pulses
time.pulse{3}.fwhm      = t0/5;          % FWHM (/length) of pulses
time.pulse{4}.fwhm      = 0;             % FWHM (/length) of pulses

