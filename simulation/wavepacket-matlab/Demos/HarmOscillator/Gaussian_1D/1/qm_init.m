% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space time

prt.disp ( '***************************************' )
prt.disp ( 'Squeezed state of a harmonic oscillator' )
prt.disp ( '***************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % Use Fourier grid
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -10;                % Lower bound of the grid
space.dof{1}.x_max = 10;                 % Upper bound of the grid
space.dof{1}.mass  = 1;                  % Mass for the kinetic energy

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 020;               % Index of final time step
time.steps.m_delta  = pi/10;             % Size of time steps 
time.steps.s_number =  100;              % Number of sub steps per time step

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width = sqrt(2);             % Width 
time.dof{1}.pos_0 = -3.0;                % Center in position representation
time.dof{1}.mom_0 =  0.0;                % Center in momentum representation

% Hamiltonian operator 
hamilt.truncate.e_min  =  0.0;           % Lower truncation of energy
hamilt.truncate.e_max  = 50.0;           % Upper truncation of energy

hamilt.pot{1,1}        = pot.taylor;     % Taylor expansion
hamilt.pot{1,1}.coeffs = [0;1];          % Force constant

% Plot time evolution of density
plots.density            = vis.contour;
plots.density.marginals  = true;
plots.density.energy     = true;
% plots.density.expect     = true;
plots.density.cnt_nlev   = [15 25];
plots.density.cnt_levels = false; 
plots.density.cnt_min    = 0.010;
plots.density.cnt_max    = 0.050;
plots.density.images     = false;

% Plot expectation values
plots.expect          = vis.expect;
plots.expect.e_min    = -2;              % Manual setting of range for energies
plots.expect.e_max    = 15;              % Manual setting of range for energies
plots.expect.p_min    = 0.99;            % Manual setting of range for population
plots.expect.p_max    = 1.01;            % Manual setting of range for population
plots.expect.errorbar = false;           % Error bars for pot./kin. energy
plots.expect.legends  = true;            % Manual setting of range for population
plots.expect.export   = true;            % Export as JPG file

% Plot spectrum
plots.spectrum = vis.spectrum;          % Power spectrum
plots.spectrum.export = true;            % Export as JPG file

