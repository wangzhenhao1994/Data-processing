% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global hamilt plots space state time

prt.disp ( '****************************************************************' )
prt.disp ( 'Dual crossing example, k=30 (Tully 1990)' )
prt.disp ( '****************************************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs      = 2;
hamilt.coupling.represent  = 'adi';
hamilt.coupling.ini_rep    = 'adi';
hamilt.coupling.ini_coeffs = [1 0];      % Initially only lower adiabatic state 

% Spatial discretization
space.dof{1} = dof.fft;                  % Using FFT grid
space.dof{1}.mass = 2000;                % Mass of paarticle
space.dof{1}.n_pts = 512;                % Number of grid points
space.dof{1}.x_min = -12.0;              % Lower bound of grid 
space.dof{1}.x_max =  12.0;              % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 040;               % Index of final time step
time.steps.m_delta  = 025;               % Size of time steps 
time.steps.s_number = 050;               % Number of sub steps per time step

% Surface hopping details (if applicable)
if isa (state,'sht.generic')
    state.rescale = 1;                   % Rescale momenta after hopping
    state.sca_nac = 0;                   % Rescale along NAC coupling vectors
end

% Hamiltonian operator
hamilt.truncate.e_min = -0.10;           % Lower truncation of energy
hamilt.truncate.e_max =  1.15;           % Upper truncation of energy

hamilt.pot{1,1} = pot.tully2;            % Double crossing model
hamilt.pot{2,2} = pot.tully2;            % Double crossing model
hamilt.pot{1,2} = pot.tully2;            % Double crossing model

% Absorbing boundary conditions
hamilt.nip{1} = nip.power;               % Negative imaginary potential
hamilt.nip{1}.exp  = 4;                  % Exponent
hamilt.nip{1}.min = -10;                 % Beginning of inner grid region
hamilt.nip{1}.max = +10;                 % End of inner grid region

hamilt.nip{2} = nip.power;               % Negative imaginary potential
hamilt.nip{2}.exp  = 4;                  % Exponent
hamilt.nip{2}.min = -10;                 % Beginning of inner grid region
hamilt.nip{2}.max = +10;                 % End of inner grid region

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width =  0.75;               % Width 
time.dof{1}.pos_0 = -7.0;                % Center in position representation
time.dof{1}.mom_0 = 30.0;                % Center in momentum representation

% Plot densities
plots.density        = vis.contour;
plots.density.expect = false;
plots.density.range  = true;             % manual setting of plotting range
plots.density.x_min = -9;
plots.density.x_max = 9;
plots.density.y_min = 10;
plots.density.y_max = 35;
plots.density.pot_max = 0.40;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max =  0.3;

