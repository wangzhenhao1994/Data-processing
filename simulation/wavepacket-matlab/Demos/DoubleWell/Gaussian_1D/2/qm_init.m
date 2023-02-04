% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it 
% into other works without any restrictions.

function qm_init
global hamilt plots space time

prt.disp ( '***********************************************' )
prt.disp ( 'Razavy symmetric double well potential' )
prt.disp ( 'Starting with Gaussian in the left well' )
prt.disp ( 'Initial energy far below threshold' )
prt.disp ( '***********************************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % using FFT grid
space.dof{1}.mass = 1/2;                 % Particle mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -7.0;               % Lower bound of grid 
space.dof{1}.x_max =  7.0;               % Upper bound of grid

% Hamiltonian operator 
hamilt.truncate.e_min  =  -15.0;         % Lower truncation of energy
hamilt.truncate.e_max  =  150.0;         % Upper truncation of energy

% Razavy Potential: beta=0.1, kappa=-7
hamilt.pot{1,1}          = pot.razavy;   % Hyperbolic potential
hamilt.pot{1,1}.modified = true;         % Use modified version
hamilt.pot{1,1}.eta = -0.7;              % prefactor of cosh
hamilt.pot{1,1}.zeta = 0.01;             % prefactor of cosh^2

% Temporal discretization and propagator
time.steps.m_start = 000;                % Index of initial time step
time.steps.m_stop  = 050;                % Index of final time step
time.steps.m_delta  = 0.2;               % Size of time steps 
time.steps.s_number   =  100;            % Number of sub steps per time step

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width = 0.34;                % Width of Gaussian
time.dof{1}.pos_0 = -4.2483;             % Centered at left minimum (well)
time.dof{1}.mom_0 = +2.8;                % Initial momentum

% Plot of densities
plots.density          = vis.contour;
plots.density.pot_max  = -15;
plots.density.pot_max  = +15;
plots.density.expect   = false;
plots.density.cnt_nlev = [30 30];

% Plot of expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 20;