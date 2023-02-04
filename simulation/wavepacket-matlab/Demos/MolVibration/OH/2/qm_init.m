% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init 
global atomic control hamilt plots space state time

state.save_export = true; % Save wavefunctions etc

prt.disp ( '***************************************************************' )
prt.disp ( 'Vibrationally state selective excitation of a Morse ' )
prt.disp ( 'oscillator resembling one of the OH bonds in water  ' )
prt.disp ( 'Laser pulse optimized for 5-photon transition: 0->5 ' )
prt.disp ( '                                                    ' )
prt.disp ( 'see M.V.Korolkov, G.K.Paramonov, and B. Schmidt     ' )
prt.disp ( 'J. Chem. Phys. 105(5), 1862-1879  (1996), Fig. 2a   ' )
prt.disp ( '***************************************************************' )

% In the paper we used a reduced mass of 1728.539, but here 1728.257
m_O = atomic.mass.O16;
m_H = atomic.mass.H1;

% Spatial discretization
space.dof{1} = dof.fft;                  % Using fft grid
space.dof{1}.mass  = m_O*m_H/(m_O+m_H);  % Reduced mass 
space.dof{1}.n_pts = 512;                % Number of grid points
space.dof{1}.x_min = 0.7;                % Lower bound of grid 
space.dof{1}.x_max = 12.0;               % Upper bound of grid

% Hamiltonian operator 
hamilt.truncate.e_min = -0.1;            % Lower truncation of energy
hamilt.truncate.e_max =  0.5;            % Upper truncation of energy

hamilt.pot{1,1}      = pot.morse;        % Harmonic oscillator
hamilt.pot{1,1}.d_e  = 0.1994;           % Dissociation energy
hamilt.pot{1,1}.r_e  = 1.821;            % Equilibrium length
hamilt.pot{1,1}.alf  = 1.189;            % Range parameter

hamilt.dip{1}{1,1} = dip.mecke;          % Mecke dipole function
hamilt.dip{1}{1,1}.r_0 = 0.6/atomic.l.A; % Length parameter: 0.6 A
hamilt.dip{1}{1,1}.q_0 = 7.85/atomic.d.D*atomic.l.A; % Charge parameter: 7.85 D/A

% Temporal discretization
time.steps.m_delta  = 1/atomic.t.fs;     % Size of time steps: 10 fs 
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 1000;              % Index of final time step
time.steps.s_number = 10;                % Number of sub steps per time step

% Calculate and save (bound) eigen states (==> qm_bound)
hamilt.eigen.stop = 00;                  % Lower index
hamilt.eigen.stop = 21;                  % Upper index

% Electric field pulses
time.pulse{1}       = efi.sin_2;         % Shape of envelope
time.pulse{1}.delay = 500/atomic.t.fs;   % Time delay of pulse center
time.pulse{1}.fwhm  = 500/atomic.t.fs;   % Full width at half maximum
time.pulse{1}.ampli = 328.5/atomic.F.MV_cm;% From MV/cm to E_H/(e*a_0)
time.pulse{1}.frequ = 3424.19/atomic.w.cm_1;% From cm-1 to E_H

% Define populations as observables
control.observe.types = 'prj';           % types of observables
control.observe.choices = {[0] [1] [2] [3] [4] [5] [6] [7:21]};
control.observe.labels  = {'|0>' '|1>' '|2>' '|3>' '|4>' '|5>' '|6>' '|7:21>'};
control.observe.targets = 1:8;           % evaluate targets at terminal time

% Initial state (==> qm_abncd)
time.ket.choice = 'pure';                % starting from a pure state
time.ket.pure = 0;                       % vibrational ground state

% Plot densities: contours or bar
if strcmpi(class(state),'wave')
    plots.density = vis.curve;           % Curve plot of wave function
elseif strcmpi(class(state),'ket')
    plots.density = vis.bar;             % Bar graph of "ket" vector
end

% Plot expectation values
if strcmpi(class(state),'wave')
    plots.expect = vis.expect;           % "normal" representation
    plots.expect.e_min = -0.05;          % Range of energy plot
    plots.expect.e_max = +0.25;          % Range of energy plot
elseif strcmpi(class(state),'ket')
    plots.expect = vis.in_out;           % input (fields u(t)) and output (observables y(t))
    plots.expect.y_min = -0.10;          % Range of population plot
    plots.expect.y_max = +1.10;          % Range of population plot
end