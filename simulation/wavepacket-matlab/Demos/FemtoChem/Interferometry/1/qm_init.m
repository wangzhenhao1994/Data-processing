% Copyright (C) 2017 - .... Burkhard Schmidt 
%               2009 - 2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global atomic hamilt plots space time

prt.disp ('**************************************')
prt.disp ('Relaxation of an wavefunction in the  ')
prt.disp ('electronic groundstate of Na2         ')
prt.disp ('**************************************')

hamilt.coupling.labels = {'X^1\Sigma_g^+'};

% Grid definition
space.dof{1}       = dof.fft;            % using FFT grid
space.dof{1}.n_pts = 256;                % number of points
space.dof{1}.x_min = 4;                  % lower boundary of the grid
space.dof{1}.x_max = 14;                 % upper boundary of the grid
space.dof{1}.mass  = atomic.mass.Na23/2; % reduced mass of 23Na2

% Temporal discretisation
time.steps.m_start  = 0;                 % index of first time step
time.steps.m_stop   = 10;                % index of last time step
time.steps.m_delta  = 10/atomic.t.fs;    % 10 fs per time step
time.steps.s_number = 1;                 % propagation time step not relevant (Chebychev)

% Hamiltonian operator
hamilt.truncate.e_min = -0.03;           % lower truncation of energy
hamilt.truncate.e_max = 0.05;            % upper truncation of energy

hamilt.pot{1,1}          = pot.interp;   % interpolate tabulated potential
hamilt.pot{1,1}.pos_conv = atomic.l.A;   % conversion factor for coordinates
hamilt.pot{1,1}.pot_conv = atomic.E.eV;  % conversion factor for energies

% Initial wave function; Gaussian around the equilibrium
time.dof{1}        = init.gauss;         % Gaussian initial wave function
time.dof{1}.width  = 0.5;                % width of Gaussian
time.dof{1}.pos_0  = 6;                  % center in position representation
time.dof{1}.mom_0  = 0;                  % center in momentum representation

% Plot densities
plots.density           = vis.contour;   % Draw contour lines
plots.density.range     = true;          % manually set ranges
plots.density.x_min     = 4;
plots.density.x_max     = 8;
plots.density.y_min     = -10;
plots.density.y_max     = +10;
plots.density.scale_dvr = 2;             % manual height of densities
plots.density.expect    = false;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 4e-4;

