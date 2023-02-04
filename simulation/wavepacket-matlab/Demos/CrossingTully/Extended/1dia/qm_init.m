% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global hamilt plots space time

prt.disp ( '****************************************************************' )
prt.disp ( 'Extended coupling example, k=10 (Tully 1990)' )
prt.disp ( '****************************************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs      = 2;          % Two coupled channels
hamilt.coupling.represent  = 'dia';      % Adiabatic representation
hamilt.coupling.ini_rep    = 'dia';      % Adiabatic initial state
hamilt.coupling.ini_coeffs = [1 0];      % Initially only lower adiabatic state populated

% Spatial discretization
space.dof{1} = dof.fft;                  % using FFT grid methods
space.dof{1}.mass = 2000;                % Mass of particle
space.dof{1}.n_pts = 768;                % Number of grid points
space.dof{1}.x_min = -15.0;              % Lower bound of grid 
space.dof{1}.x_max =  10.0;              % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 040;               % Index of final time step
time.steps.m_delta  = 100.0;             % Size of time steps 
time.steps.s_number = 050;               % Number of sub steps per time step

% Hamiltonian operator 
hamilt.truncate.e_min    = -0.20;        % Lower truncation of energy
hamilt.truncate.e_max    =  0.50;        % Upper truncation of energy

hamilt.pot{1,1}   = pot.tully3;          % Extended Crossing
hamilt.pot{2,2}   = pot.tully3;          % Extended Crossing
hamilt.pot{1,2}   = pot.tully3;          % Extended Crossing

% Absorbing boundary conditions
hamilt.nip{1} = nip.power;               % Negative imaginary potential
hamilt.nip{1}.exp  = 4;                  % Exponent
hamilt.nip{1}.min = -13;                 % Beginning of inner grid region
hamilt.nip{1}.max = +08;                 % End of inner grid region

hamilt.nip{2} = nip.power;               % Negative imaginary potential
hamilt.nip{2}.exp  = 4;                  % Exponent
hamilt.nip{2}.min = -13;                 % Beginning of inner grid region
hamilt.nip{2}.max = +08;                 % End of inner grid region

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width =   0.75;              % Width 
time.dof{1}.pos_0 = -10.0;               % Center in position representation
time.dof{1}.mom_0 =  10.0;               % Center in momentum representation

% Plot densities
plots.density        = vis.curve;
plots.density.pot_min = -0.2; 
plots.density.pot_max = +0.4;

% plots.density        = vis.contour;
% plots.density.expect = false;
% plots.density.range  = true;             % manual setting of plotting range
% plots.density.x_min = -14;
% plots.density.x_max = 9;
% plots.density.y_min = -20;
% plots.density.y_max = 40;
% plots.density.pot_min = -0.2; 
% plots.density.pot_max = +0.2; 

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 0.3;