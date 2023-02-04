% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space time

prt.disp ( '***************************************************************' )
prt.disp ( 'Coherent state of a harmonic oscillator' )
prt.disp ( '***************************************************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % Using FFT grid method
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max =  10.0;              % Upper bound of grid
space.dof{1}.mass = 1;                   % Mass

% Temporal discretization
time.steps.m_start = 000;                % Index of initial time step
time.steps.m_stop  = 020;                % Index of final time step
time.steps.m_delta = pi/10;              % Size of time steps 
time.steps.s_number = 100;               % Number of sub steps per time step

% Initial state
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width = sqrt(1/2);           % Width 
time.dof{1}.pos_0 = -5.0;                % Center in position representation
time.dof{1}.mom_0 =  0.0;                % Center in momentum representation

% Truncate Hamiltonian operator 
hamilt.truncate.e_min  =  0.0;           % Lower truncation of energy
hamilt.truncate.e_max  = 75.0;           % Upper truncation of energy

% Potential energy
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
plots.spectrum = vis.spectrum;           % Power spectrum
plots.spectrum.export = true;            % Export as JPG file

% plots.density = vis.curve;               % Contour plot of the Wigner function
% plots.density.represent = 'dvr';         % Position or momentum space
% plots.density.complex = 'abs2';          % Real or Imag or Abs2
% plots.density.energy = true;             % Show also energy function
% plots.density.range    = true;           % Manual setting of plotting range
% plots.density.x_min = -8;
% plots.density.x_max = 8;
% plots.density.y_min = -8;
% plots.density.y_max = 8;

% plots.density = vis.flux;

% plots.density = vis.surface;
% plots.density.srf_view  = [75 60];       % az el
% plots.density.srf_light = [25 75];       % az el
% plots.density.col_map = 'default';        
% plots.density.energy = true;
% plots.density.marginals = true;
% plots.density.range    = true;           % Manual setting of plotting range
% plots.density.x_min = -8;
% plots.density.x_max = 8;
% plots.density.y_min = -8;
% plots.density.y_max = 8;
% plots.density.z_min = 0;
% plots.density.z_max = 50;


