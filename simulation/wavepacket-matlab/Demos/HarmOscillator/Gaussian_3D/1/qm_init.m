% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space time

prt.disp ( '************************************************' )
prt.disp ( 'Squeezed/coherent state of a harmonic oscillator' )
prt.disp ( '************************************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % Using FFT grid 
space.dof{1}.n_pts = 032;                % Number of grid points
space.dof{1}.x_min = 05.0;               % Lower bound of grid 
space.dof{1}.x_max = 15.0;               % Upper bound of grid
space.dof{1}.mass = 1.0;                 % Mass

space.dof{2} = dof.fft;                  % Using FFT grid 
space.dof{2}.n_pts = 064;                % Number of grid points
space.dof{2}.x_min = 15.0;               % Lower bound of grid 
space.dof{2}.x_max = 25.0;               % Upper bound of grid
space.dof{2}.mass = 1.0;                 % Mass

space.dof{3} = dof.fft;                  % Using FFT grid 
space.dof{3}.n_pts = 064;                % Number of grid points
space.dof{3}.x_min = 25.0;               % Lower bound of grid 
space.dof{3}.x_max = 35.0;               % Upper bound of grid
space.dof{3}.mass = 1.0;                 % Mass

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 050;               % Index of final time step
time.steps.m_delta  = pi/25;             % Size of time steps 
time.steps.s_number = 0100;              % Number of sub steps per time step

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width = 1.5*sqrt(1/2);       % Width: squeezed state
time.dof{1}.pos_0 = 10.0;                % Center in position representation
time.dof{1}.mom_0 =  0.0;                % Center in momentum representation

time.dof{2}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{2}.width = sqrt(1/2);           % Width: coherent state
time.dof{2}.pos_0 =  22.5;               % Center in position representation
time.dof{2}.mom_0 =  +0.0;               % Center in momentum representation

time.dof{3}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{3}.width = sqrt(1/2);           % Width: coherent state
time.dof{3}.pos_0 =  30.0;               % Center in position representation
time.dof{3}.mom_0 =  -2.5;               % Center in momentum representation

% Hamiltonian operator 
hamilt.truncate.e_min  =   0.0;          % Lower truncation of energy
hamilt.truncate.e_max  = 100.0;          % Upper truncation of energy

hamilt.pot{1,1} = pot.taylor;            % Taylor series: Harmonic oscillator
hamilt.pot{1,1}.coeffs = [0 0 0;1 1 1];  % Force constants
hamilt.pot{1,1}.hshift = [10  20  30 ];  % Minimum positions

hamilt.nip{1}     = nip.power;           % Negativ imaginary potential
hamilt.nip{1}.exp = [2  2  2 ];          % Exponent
hamilt.nip{1}.min = [02 12 22];          % Lower bound
hamilt.nip{1}.max = [18 28 38];          % Upper bound

% Plot density using isosurfaces
plots.density           = vis.surface;   % 3D surface plot
plots.density.represent = 'dvr';
plots.density.srf_view  = [60 45];       % View point for surface plot: [az el]

% plots.density           = vis.reduced_1d;  % 3 Wigner contour plots
% plots.density.represent = 'wig';

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_min = 00;                 % Manually set range for energies
plots.expect.e_max = 10;                 % Manually set range for energies
