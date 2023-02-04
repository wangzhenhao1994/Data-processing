% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space time

prt.disp ( '*******************************************************************************' )
prt.disp ( 'Free particle with momentum: Wavepacket translation and dispersion' )
prt.disp ( '*******************************************************************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % using fft grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max = +20.0;              % Upper bound of grid

% Spatial discretization
space.dof{2} = dof.fft;                  % using fft grid
space.dof{2}.mass = 1;                   % mass
space.dof{2}.n_pts = 128;                % Number of grid points
space.dof{2}.x_min = -10.0;              % Lower bound of grid 
space.dof{2}.x_max = +20.0;              % Upper bound of grid

% Temporal discretization and propagator
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 063;               % Index of final time step
time.steps.m_delta  = 0.1;               % Size of time steps 
time.steps.s_number = 020;               % Number of sub steps per time step

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width = sqrt(1/2) * 1;       % Width 
time.dof{1}.pos_0 = -5.0;                % Center in position representation
time.dof{1}.mom_0 = +2.0;                % Center in momentum representation

time.dof{2}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{2}.width = sqrt(1/2) * 1;       % Width 
time.dof{2}.pos_0 = -5.0;                % Center in position representation
time.dof{2}.mom_0 = +2.0;                % Center in momentum representation

% Hamiltonian operator 
hamilt.truncate.e_min  =   0.0;          % Lower truncation of energy
hamilt.truncate.e_max  =  30.0;          % Upper truncation of energy

% Absorbing boundary conditions
hamilt.nip{1}     = nip.power;           % Negative imaginary potential
hamilt.nip{1}.exp = [2 2];               % Exponent
hamilt.nip{1}.min = [-07 -07];           % Beginning of inner grid region
hamilt.nip{1}.max = [+15 +15];           % End of inner grid region

% Plot densities
plots.density           = vis.surface;   % make a neat 3D plot
plots.density.energy    = false;         % Densities with phases (color-coded)
plots.density.srf_view  = [45 15];       % Specify az/el viewing angles
plots.density.represent = 'dvr';

% Plot expectation values
plots.expect      = vis.expect;
plots.expect.e_max=5.0;
