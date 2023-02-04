% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global atomic hamilt plots space time

prt.disp ( '*****************************************************************' )
prt.disp ( 'Retinal isomerization: model in 1 dimension' )
prt.disp ( 'S.Hahn, G.Stock, J. Phys. Chem. B 104(6),1149 (2000)' )
prt.disp ( '*****************************************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs      = 2;
hamilt.coupling.represent  = 'adi';
hamilt.coupling.ini_rep    = 'adi';
hamilt.coupling.ini_coeffs = [0 1];      % Initially only upper adiabatic state populated

% Spatial discretization
space.dof{1}       = dof.fft;            % using FFT grid
space.dof{1}.mass  = 56198.347;          % effective mass
space.dof{1}.n_pts = 384;                % Number of grid points
space.dof{1}.x_min = -pi;                % Lower bound of grid 
space.dof{1}.x_max =  pi;                % Upper bound of grid

space.dof{2}       = dof.fft;            % using FFT grid
space.dof{2}.mass  = 143.16;             % effective mass
space.dof{2}.n_pts = 64;                 % Number of grid points
space.dof{2}.x_min = -4;                 % Lower bound of grid 
space.dof{2}.x_max =  4;                 % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 050;               % Index of final time step
time.steps.m_delta  = 100.0;             % Size of time steps 
time.steps.s_number = 100;               % Number of sub steps per time step

% Hamiltonian operator 
hamilt.truncate.e_min    =  -0.1;        % Lower truncation of energy
hamilt.truncate.e_max    =  +0.5;        % Upper truncation of energy

hamilt.pot{1,1}         = pot.retinal;         % Hahn&Stock: Retinal potential
hamilt.pot{1,1}.shift   = 0.0  / atomic.E.eV;  % Vertical offset
hamilt.pot{1,1}.barrier = 3.6  / atomic.E.eV;  % Rotational barrier
hamilt.pot{1,1}.omega   = 0.19 / atomic.E.eV;  % Vibrational frequency
hamilt.pot{1,1}.kappa   = 0.0  / atomic.E.eV;  % Linear parameter

hamilt.pot{2,2}         = pot.retinal;         % Hahn&Stock: Retinal potential
hamilt.pot{2,2}.shift   =  2.48 / atomic.E.eV; % Vertical offset
hamilt.pot{2,2}.barrier = -1.09 / atomic.E.eV; % Rotational barrier
hamilt.pot{2,2}.omega   =  0.19 / atomic.E.eV; % Vibrational frequency
hamilt.pot{2,2}.kappa   = -0.10 / atomic.E.eV; % Linear parameter

hamilt.pot{1,2}         = pot.retinal;         % Hahn&Stock: Retinal potential
hamilt.pot{1,2}.lambda  = 0.19 / atomic.E.eV;  % Linear coupling parameter

% Initial wave function
time.dof{1}        = init.gauss;         % Gaussian-shaped wavepacket
time.dof{1}.width  = 0.25;               % Width 
time.dof{1}.pos_0  = 0;                  % Center in position representation
time.dof{1}.mom_0  = 0;                  % Center in momentum representation

time.dof{2}        = init.gauss;         % Gaussian-shaped wavepacket
time.dof{2}.width  = 1;                  % Width 
time.dof{2}.pos_0  = 0;                  % Center in position representation
time.dof{2}.mom_0  = 0;                  % Center in momentum representation

% Plot densities
plots.density          = vis.surface;
plots.density.srf_view = [32 32];        % View point for surface plot: [az el]
plots.density.energy   = true;

% Plot expectation values
plots.expect        = vis.expect;
plots.expect.e_max  = 0.15;
