% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global hamilt plots space time

prt.disp ( '******************************************' )
prt.disp ( 'Generic E x e conical intersection with   ' )
prt.disp ( 'linear and quadratic Jahn-Teller coupling:' )
prt.disp ( 'Dynamics on warped, lower state only' )
prt.disp ( 'Localized wavepacket without momentum' )
prt.disp ( '******************************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % using FFT grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 064;                % Number of grid points
space.dof{1}.x_min = -7.5;               % Lower bound of grid 
space.dof{1}.x_max =  7.5;               % Upper bound of grid

% Spatial discretization
space.dof{2} = dof.fft;                  % using FFT grid
space.dof{2}.mass = 1;                   % mass
space.dof{2}.n_pts = 064;                % Number of grid points
space.dof{2}.x_min = -7.5;               % Lower bound of grid 
space.dof{2}.x_max =  7.5;               % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 050;               % Index of final time step
time.steps.m_delta  = 0.05;              % Size of time steps 
time.steps.s_number = 0050;              % Number of sub steps per time step

% Hamiltonian operator 
hamilt.truncate.e_min = -15.0;           % Lower truncation of dia. energy
hamilt.truncate.e_max = +25.0;           % Upper truncation of dia. energy

hamilt.pot{1,1}       = pot.con_int;     % Conical intersection
hamilt.pot{1,1}.omega =  5.0;            % Harmonic frequency
hamilt.pot{1,1}.kappa = 10.0;            % Linear JT coupling
hamilt.pot{1,1}.gamma =  1.0;            % Quadratic JT coupling

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width =  0.5;                % Width 
time.dof{1}.pos_0 =  2.5;                % Center in position representation
time.dof{1}.mom_0 =  0.0;                % Center in momentum representation

time.dof{2}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{2}.width =  0.5;                % Width 
time.dof{2}.pos_0 =  0.0;                % Center in position representation
time.dof{2}.mom_0 =  0.0;                % Center in momentum representation

% Plot densities
plots.density        = vis.contour;
plots.density.expect = false;
plots.density.range  = true;             % Manual setting of plotting range
plots.density.x_min  = -5;
plots.density.x_max  =  5;
plots.density.y_min  = -5;
plots.density.y_max  =  5;

% plots.density           = vis.surface;
% plots.density.srf_view  =[45 75];   % Perspective view 

% Plot expectation values
plots.expect            = vis.expect;
plots.expect.e_max      = 4;             % manual setting of upper energy bond