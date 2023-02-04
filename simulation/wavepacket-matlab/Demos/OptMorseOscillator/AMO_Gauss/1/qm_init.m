% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic control hamilt plots space state time 

prt.disp ( '*****************************************************************' )
prt.disp ( 'Optimal control of a bond length for a model Morse  ' )
prt.disp ( 'oscillator resembling one of the OH bonds in water  ' )
prt.disp ( '                                                    ' )
prt.disp ( 'see Figs. 1-4 of: W. Zhu and H. Rabitz              ' )
prt.disp ( 'J. Chem. Phys. 109(2), 385-391  (1998)              ' )
prt.disp ( 'http://dx.doi.org/10.1063/1.476575                  ' )
prt.disp ( '*****************************************************************' )

% Save wavefunctions etc  to file(s)
state.save_export = true; 

% Isotopic masses
m_O = atomic.mass.O16;
m_H = atomic.mass.H1;

% Spatial discretization ==> qm_bound
space.dof{1}       = dof.fft;            % Using fft grid
space.dof{1}.mass  = m_O*m_H/(m_O+m_H);  % Reduced mass (OH radical)
space.dof{1}.n_pts = 512;                % Number of grid points
space.dof{1}.x_min =  0.7;               % Lower bound of grid 
space.dof{1}.x_max = 10.0;               % Upper bound of grid     

hamilt.amo{1}       = amo.gauss;         % Choose Gaussian-shaped AMO
hamilt.amo{1}.pos_0 = 2.5;               % Center of Gaussian
hamilt.amo{1}.width = 1/(2*25);          % Width parameter of Gaussian
hamilt.amo{1}.label = 'Gaussian|2.5|25';  % Labeling this AMO

% Hamiltonian operator 
hamilt.truncate.e_min = -0.1;            % Lower truncation of energy
hamilt.truncate.e_max = +1.0;            % Upper truncation of energy

hamilt.pot{1,1}      = pot.morse;        % Morse oscillator
hamilt.pot{1,1}.d_e  = 0.1994;           % Dissociation energy
hamilt.pot{1,1}.r_e  = 1.821;            % Equilibrium length
hamilt.pot{1,1}.alf  = 1.189;            % Range parameter

hamilt.dip{1}{1,1}     = dip.mecke;      % Mecke dipole function
hamilt.dip{1}{1,1}.r_0 = 0.6/atomic.l.A; % Length parameter: 0.6 A
hamilt.dip{1}{1,1}.q_0 = 7.85/atomic.d.D*atomic.l.A; % Charge parameter: 7.85 D/A

% Temporal discretization
time.steps.m_delta  = 10/atomic.t.fs;    % Size of time steps: 10 fs 
time.steps.m_start  = 0000;              % Index of initial time step
time.steps.m_stop   = 0317;              % Index of final time step
time.steps.s_number = 300;               % Number of sub steps per time step

% Calculate and save (bound) eigen states (==> qm_bound)
hamilt.eigen.stop = 00;                  % Lower index
hamilt.eigen.stop = 21;                  % Upper index

% Electric field pulses: Initial guess
time.pulse{1}        = efi.recta;        % Shape of envelope
time.pulse{1}.delay  = 1585/atomic.t.fs; % Time delay of pulse center
time.pulse{1}.fwhm   = 3170/atomic.t.fs; % Pulse length
time.pulse{1}.ampli  = 0.02;             % Amplitude of electric field         

time.frog.choice = 'none';               % Frequency resolved optical gating
time.frog.zoom = 15;                     % Zoom factor for frequency axes 

% Initial state (==> qm_abncd)
time.ket.choice = 'pure';                % starting from a pure state
time.ket.pure = 0;                       % vibrational ground state

% Define outputs: projections as observables
control.observe.types = 'amo';           % types of observables
control.observe.choices = {1};
control.observe.targets = 1;             % choose control targets 

% Optimal control 
control.optimal.terminal = 1;            % Which observable to be optimized
control.optimal.tolerance = 1e-20;       % Threshold terminating iteration
control.optimal.max_iter = 040;          % Max. number of iteration steps
control.optimal.alpha = 1.00;            % Penalty factor for laser fluence
control.optimal.eta  = 1.00;             % Calculate optimal backward field
control.optimal.zeta = 1.00;             % Calculate optimal forward field
control.optimal.order = 2;               % Error order for optimal fields
control.optimal.prefactor = 'current';   % How to calculate prefactor (C2)
control.optimal.fb_test = false;         % Only testing forward/backward

% Plotting options
control.plot.uxy = true;                 % Plot u(t), x(t), y(t)
control.plot.j12 = true;                 % Plot j_1(t), j_2(t), and total
control.plot.psd = true;                 % Plot power spectral density
control.plot.mov = true;                 % Animation of u(t), x(t), y(t)

% Plot u(t), x(t), y(t)
plots.control = vis.uxy;


