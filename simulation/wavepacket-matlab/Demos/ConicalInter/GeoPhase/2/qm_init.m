% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global hamilt plots space time

prt.disp ( '******************************************' )
prt.disp ( 'Generic E x e conical intersection example' )
prt.disp ( 'with linear Jahn-Teller coupling included:' )
prt.disp ( 'Dynamics on lower (Mexican hat) state only' )
prt.disp ( 'Neglecting geometric (Berry) phase effect!' )
prt.disp ( '******************************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % using fft grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = -4.5;               % Lower bound of grid 
space.dof{1}.x_max =  4.5;               % Upper bound of grid

% Spatial discretization
space.dof{2} = dof.fft;                  % using fft grid
space.dof{2}.mass = 1;                   % mass
space.dof{2}.n_pts = 64;                 % Number of grid points
space.dof{2}.x_min = -4.5;               % Lower bound of grid 
space.dof{2}.x_max =  4.5;               % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 100;               % Index of final time step
time.steps.m_delta  = 0.05;              % Size of time steps 
time.steps.s_number = 50;                % Number of sub steps per time step

% Hamiltonian operator 
hamilt.truncate.e_min    = -200.0;       % Lower truncation of dia. energy
hamilt.truncate.e_max    = +200.0;       % Upper truncation of dia. energy

hamilt.pot{1,1}       = pot.con_int;     % Conical intersection
hamilt.pot{1,1}.omega =  5.0;            % Harmonic frequency
hamilt.pot{1,1}.kappa = 10.0;            % Linear JT coupling
hamilt.pot{1,1}.gamma =  0.0;            % Quadratic JT coupling

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width =  0.5;                % Width 
time.dof{1}.pos_0 = -2.0;                % Center in position representation
time.dof{1}.mom_0 =  0.0;                % Center in momentum representation

time.dof{2}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{2}.width =  0.5;                % Width 
time.dof{2}.pos_0 =  0.0;                % Center in position representation
time.dof{2}.mom_0 =  0.0;                % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density           = vis.surface;
plots.density.srf_view  = [30 60];       % Side view 
plots.density.energy    = true;          % With energy surfaces, densities color-coded
% plots.density.surface.view  = [90 -90];  % Bottom view
% plots.density.srf_look  = [ 1  1];       % Shading and lighting
% plots.density.srf_lite  = [90 -90];      % Angle of light for surface plot: [az el]
% plots.density.energy     = false;        % No energy surfaces, only densities

% plots.density         = vis.contour;
% plots.density.eexpect = false;
% plots.density.range   = true;  % manual setting of ranges of the plots
% plots.density.x_min   = -3.5;
% plots.density.x_max   = +3.5;
% plots.density.y_min   = -3.5;
% plots.density.y_max   = +3.5;

plots.expect = vis.expect;
plots.expect.e_min   = -15;              % manual setting of minimum of energy expct. plot
plots.expect.e_max   = +30;              % manual setting of maximum of energy expct. plot
