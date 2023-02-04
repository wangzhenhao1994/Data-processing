% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots space time

prt.disp ( '************************************************' )
prt.disp ( 'Fulvene torsionial dynamics (electronic B state)' )
prt.disp ( 'See http://dx.doi.org/10.1002/cphc.200600543    ' )
prt.disp ( '************************************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % using fft grid
space.dof{1}.mass = 1/(2.44/atomic.E.meV);  % Reduced moment of inertia: 1/2.44 meV
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -pi;                % Lower bound of grid 
space.dof{1}.x_max =  pi;                % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 075;               % Index of final time step
time.steps.m_delta  = 5/atomic.t.fs;     % Size of time steps: 05 fs
time.steps.s_number =  50;               % Number of sub steps per time step

% Hamiltonian operator 
hamilt.truncate.e_min    = -0.015;       % Lower truncation of energy
hamilt.truncate.e_max    = +0.050;       % Upper truncation of energy

hamilt.pot{1,1}      = pot.pendulum;     % Intramolecular torsion
hamilt.pot{1,1}.zeta = 106*2.44/atomic.E.meV; % Barrier height: 106*2.44 meV

% Initial wave function
time.dof{1}          = init.pendulum1;   % Eigenstate of plane pendulum
time.dof{1}.parity   = 'c';              % Parity: cosine=even
time.dof{1}.order    = 0;                % Order (quantum number)
time.dof{1}.multiple = 2;                % Multiplicity: Double well
time.dof{1}.barrier  = 856*2.44/atomic.E.meV*space.dof{1}.mass; % 856*2.44 meV
time.dof{1}.shift    = 0;                % Potential shift: None

% Plot densities
plots.density         = vis.contour;     % Contour plot of the eigenstates
plots.density.expect  = false;
plots.density.pot_min = -10e-3;          % adjust plot energy range
plots.density.pot_max = +04e-3;

% Plot expectation values
plots.expect       = vis.expect;
