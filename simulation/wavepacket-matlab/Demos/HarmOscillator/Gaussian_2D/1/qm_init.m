% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space time

prt.disp ( '***************************************' )
prt.disp ( 'Coherent state of a harmonic oscillator' )
prt.disp ( '***************************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % Using FFT grid 
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -15.0;              % Lower bound of grid 
space.dof{1}.x_max = +15.0;              % Upper bound of grid
space.dof{1}.mass = 1.0;                 % Mass

space.dof{2} = dof.fft;                  % Using FFT grid 
space.dof{2}.n_pts = 128;                % Number of grid points
space.dof{2}.x_min = -15.0;              % Lower bound of grid 
space.dof{2}.x_max = +15.0;              % Upper bound of grid
space.dof{2}.mass = 1.0;                 % Mass

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 040;               % Index of final time step
time.steps.m_delta  = pi/20;             % Size of time steps 
time.steps.s_number = 0100;              % Number of sub steps per time step

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width = sqrt(1/2);           % Width 
time.dof{1}.pos_0 = -5.0;                % Center in position representation
time.dof{1}.mom_0 =  0.0;                % Center in momentum representation

time.dof{2}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{2}.width = sqrt(1/2);           % Width 
time.dof{2}.pos_0 = -5.0;                % Center in position representation
time.dof{2}.mom_0 = -5.0;                % Center in momentum representation

% Hamiltonian operator 
hamilt.truncate.e_min  =   0.0;          % Lower truncation of energy
hamilt.truncate.e_max  = 100.0;          % Upper truncation of energy

hamilt.pot{1,1}        = pot.taylor;     % Taylor series: Harmonic oscillator
hamilt.pot{1,1}.coeffs = [0 0; 1 1];     % Force constants

hamilt.nip{1}      = nip.power;          % Negativ imaginary potential
hamilt.nip{1}.exp  = [2 2];              % Exponent
hamilt.nip{1}.min  = [-12.0 -12.0];      % Lower bound
hamilt.nip{1}.max  = [+12.0 +12.0];      % Upper bound

% Plot time evolution of density
plots.density           = vis.contour;   % 2D contour plot
plots.density.energy    = true;
plots.density.marginals = true;
plots.density.represent = 'dvr';

% plots.density           = vis.surface;   % 3D surface plot
% plots.density.represent = 'dvr';
% plots.density.srf_view  = [30 75];       % View point for surface plot: [az el]
% plots.density.energy    = true;

% plots.density           = vis.flux;      % 2D quiver plot

% plots.density = vis.reduced_1d;
% plots.density.represent = 'wig';

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_min = 00;                 % Manually set range for energies
plots.expect.e_max = 50;                 % Manually set range for energies
