% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global hamilt plots space time

prt.disp ( '********************************************' )
prt.disp ( 'Spin-boson system in 1 dimensions' )
prt.disp ( 'Parameters from Monika Hejjas thesis' )
prt.disp ( 'HU Berlin, Dept. of Physics, 2004' )
prt.disp ( '********************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs      = 2;
hamilt.coupling.represent  = 'adi';
hamilt.coupling.ini_rep    = 'dia';
hamilt.coupling.ini_coeffs = [1 0];      % Initially left diabatic state populated

% Spatial discretization
space.dof{1}       = dof.fft;            % using fft grid
space.dof{1}.mass  = 1;                  % mass
space.dof{1}.n_pts = 1024;               % Number of grid points
space.dof{1}.x_min = -90.0;              % Lower bound of grid 
space.dof{1}.x_max = 120.0;              % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 150;               % Index of final time step
time.steps.m_delta  = 0.5;               % Size of time steps 
time.steps.s_number = 0200;              % Number of sub steps per time step

% Hamiltonian operator 
hamilt.truncate.e_min    = -10.0;        % Lower truncation of dia. energy
hamilt.truncate.e_max    = +50.0;        % Upper truncation of dia. energy

for m = 1:hamilt.coupling.n_eqs
    for n = m:hamilt.coupling.n_eqs
        hamilt.pot{m,n} = pot.con_int;   % Conical intersection
        hamilt.pot{m,n}.omega = 00.01;   % Harmonic frequency
        hamilt.pot{m,n}.gamma = 00.50;   % spin coupling
        hamilt.pot{m,n}.kappa = sqrt(0.05); % asymmetry parameter
        hamilt.pot{m,n}.delta = 2*7.44;  % energy gap
    end
end

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width = 2.236068;            % Width 
time.dof{1}.pos_0 = -71.5;               % Center in position representation
time.dof{1}.mom_0 =  00.0;               % Center in momentum representation

% Plot densities
plots.density        = vis.contour;
plots.density.expect = false;
plots.density.energy = true;

% Plot expectation values
plots.expect            = vis.expect;
plots.expect.e_max      = 20;            % manual setting of upper energy bond