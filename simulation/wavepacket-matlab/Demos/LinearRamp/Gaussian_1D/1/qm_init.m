% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space time

prt.disp ( '*******************************************************************************' )
prt.disp ( 'Gaussian packet scattered from linear ramp' )
prt.disp ( '*******************************************************************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % using fft grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 192;                % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max =  20.0;              % Upper bound of grid

% Temporal discretization and propagator
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 020;               % Index of final time step
time.steps.m_delta  = 0.5;               % Size of time steps 
time.steps.s_number = 050;               % Number of sub steps per time step

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width = +1.0;                % Width 
time.dof{1}.pos_0 = -2.0;                % Center in position representation
time.dof{1}.mom_0 = +4.0;                % Center in momentum representation

% Hamiltonian operator 
hamilt.truncate.e_min  = -10.0;          % Lower truncation of energy
hamilt.truncate.e_max  = +30.0;          % Upper truncation of energy

hamilt.pot{1,1}        = pot.taylor;     % Taylor series: Linear potential
hamilt.pot{1,1}.coeffs = 1;              % Slope parameter

% Absorbing boundary conditions
hamilt.nip{1}      = nip.power;          % Negative imaginary potential
hamilt.nip{1}.exp  = 4;                  % Exponent
hamilt.nip{1}.min = -07;                 % Beginning of inner-grid region
hamilt.nip{1}.max =  15;                 % End of inner-grid region


% Plots of densities
plots.density       = vis.contour;       % Contour plot of Wigner transform

% Plots of expectation values
plots.expect        = vis.expect;
plots.expect.e_max  = 20;                % manually set ranges for energy plot
