% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots space time

prt.disp ( '*****************************************************************' )
prt.disp ( 'Optimal control of a bond length for a model Morse  ' )
prt.disp ( 'oscillator resembling one of the OH bonds in water  ' )
prt.disp ( '                                                    ' )
prt.disp ( 'see Figs. 1-4 of: W. Zhu and H. Rabitz              ' )
prt.disp ( 'J. Chem. Phys. 109(2), 385-391  (1998)              ' )
prt.disp ( 'http://dx.doi.org/10.1063/1.476575                  ' )
prt.disp ( '*****************************************************************' )

% Isotopic masses
m_O = atomic.mass.O16;
m_H = atomic.mass.H1;

% Spatial discretization ==> qm_bound
space.dof{1}       = dof.fft;            % Using fft grid
space.dof{1}.mass  = m_O*m_H/(m_O+m_H);  % Reduced mass (OH radical)
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = 0.7;                % Lower bound of grid 
space.dof{1}.x_max = 5.0;                % Upper bound of grid     

hamilt.amo{1}       = amo.gauss;         % Choose Gaussian-shaped AMO
hamilt.amo{1}.pos_0 = 2.5;               % Center of Gaussian
hamilt.amo{1}.width = 1/(2*25);          % Width parameter of Gaussian
hamilt.amo{1}.label = 'Gaussian|2.5|25'; % Labeling this AMO

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

% Electric field loaded from file
time.pulse{1}        = efi.inter;        % Interpolate from file
time.pulse{1}.ampli  = 1;                % Whole pulse from file
time.pulse{1}.file   = 'ket_optimal_40.dat';   % file to load the field from
time.pulse{1}.method = 'spline';         % Spline interpolation

% Initial wave function
time.dof{1}     = init.morse;            % Ground state of Morse oscillator
time.dof{1}.d_e = hamilt.pot{1,1}.d_e;   % data copied from Morse potential
time.dof{1}.r_e = hamilt.pot{1,1}.r_e;   % data copied from Morse potential
time.dof{1}.alf = hamilt.pot{1,1}.alf;   % data copied from Morse potential
time.dof{1}.m_r = space.dof{1}.mass;     % data copied from grid definition
time.dof{1}.n_q = 0;                     % ground state

% Plot densities
plots.density = vis.curve;

% Plot expectation values
plots.expect = vis.expect;