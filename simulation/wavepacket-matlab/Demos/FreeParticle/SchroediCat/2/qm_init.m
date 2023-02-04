% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space time

prt.disp ( '***********************' )
prt.disp ( 'Schroedingers cat state' )
prt.disp ( '***********************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % using fft grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -15.0;              % Lower bound of grid 
space.dof{1}.x_max = +15.0;              % Upper bound of grid

% Temporal discretization band propagator
time.steps.m_start = 000;                % Index of initial time step
time.steps.m_stop  = 010;                % Index of final time step
time.steps.m_delta  = pi/10;             % Size of time steps 
time.steps.s_number   =  100;            % Number of sub steps per time step

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.pos_0 = [  0; 0 ];           % Center in position representation
time.dof{1}.mom_0 = [ -2; 2 ];           % Center in momentum representation
time.dof{1}.width = [  1; 1 ];           % Width 

% Hamiltonian operator 
hamilt.truncate.e_min  =   0.0;          % Lower truncation of energy
hamilt.truncate.e_max  =  30.0;          % Upper truncation of energy

% Plot densities
plots.density       = vis.contour;       % Contour plot
plots.density.range = true;              % Manual setting of plotting range
plots.density.x_min = -12;
plots.density.x_max = +12;
plots.density.y_min = -04;
plots.density.y_max = +04;
plots.density.expect = false;

% Plot expectation values
plots.expect       = vis.expect;      
plots.expect.e_max = 2.5;
