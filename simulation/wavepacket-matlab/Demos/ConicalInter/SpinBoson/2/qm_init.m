% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global hamilt plots space time

prt.disp ( '********************************************' )
prt.disp ( 'Spin-boson system in 2 dimensions: adiabatic' )
prt.disp ( '********************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs      = 2;
hamilt.coupling.represent  = 'adi';
hamilt.coupling.ini_rep    = 'dia';
hamilt.coupling.ini_coeffs = [0 1];      % Initially left diabatic state populated

% Spatial discretization
space.dof{1}       = dof.fft;            % using FFT grid
space.dof{1}.mass  = 1;                  % mass
space.dof{1}.n_pts = 192;                % Number of grid points
space.dof{1}.x_min = -12.0;              % Lower bound of grid 
space.dof{1}.x_max =  12.0;              % Upper bound of grid

% Spatial discretization
space.dof{2}       = dof.fft;            % using FFT grid
space.dof{2}.mass  = 1;                  % mass
space.dof{2}.n_pts = 192;                % Number of grid points
space.dof{2}.x_min = -12.0;              % Lower bound of grid 
space.dof{2}.x_max =  12.0;              % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 042;               % Index of final time step
time.steps.m_delta  = 0.1;               % Size of time steps 
time.steps.s_number = 0150;              % Number of sub steps per time step

% Hamiltonian operator 
hamilt.truncate.e_min    = -500.0;       % Lower truncation of dia. energy
hamilt.truncate.e_max    = +500.0;       % Upper truncation of dia. energy

for m = 1:hamilt.coupling.n_eqs
    for n = m:hamilt.coupling.n_eqs
        hamilt.pot{m,n} = pot.con_int;   % Conical intersection
        hamilt.pot{m,n}.omega= 03.0;     % Harmonic frequency
        hamilt.pot{m,n}.kappa= 06.0;     % Linear JT coupling
        hamilt.pot{m,n}.gamma= 00.0;     % Quadratic JT coupling
    end
end

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width =  0.5;                % Width 
time.dof{1}.pos_0 = -5.0;                % Center in position representation
time.dof{1}.mom_0 =  0.0;                % Center in momentum representation

time.dof{2}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{2}.width =  2.5;                % Width 
time.dof{2}.pos_0 =  0.0;                % Center in position representation
time.dof{2}.mom_0 =  0.0;                % Center in momentum representation

% Plot densities
plots.density           = vis.contour;
plots.density.expect    = false;
% plots.density           = vis.surface;
% plots.density.srf_view  = [15 45];       % Perspective view (with potentials)

% Plot expectation values
plots.expect            = vis.expect;
plots.expect.e_max      = 150;           % manual setting of upper energy bond
