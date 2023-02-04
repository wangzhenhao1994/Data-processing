% Copyright (C) 2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global atomic hamilt space time

prt.disp( '***********************************************' )
prt.disp( 'Model system for testing spherical coordinate  ' )
prt.disp( 'calculations. See model I in                   ' )
prt.disp( ' J.Chem.Phys 95, 7392 (1991)                   ' )
prt.disp( '***********************************************' )

% Only one angular degree of freedom. We set R to 1 and the mass such
% that hbar^2 / 2mR^2 = 1.47822 cm^-1
space.dof{1}       = dof.legendre;       % grid for expansion in Legendre polynomials
space.dof{1}.label = 'cos \Theta';
space.dof{1}.R_0   = 1;                  % constant value for R
space.dof{1}.m_0   = 0;                  % minor quantum number fixed to 0
space.dof{1}.l_max = 50;                 % maximum angular momentum
space.dof{1}.mass  = 1 / (2 * 1.47822 / atomic.w.cm_1); % adjusted mass

% Temporal discretization
time.steps.m_start = 000;                % Index of initial time step
time.steps.m_stop  = 150;                % Index of final time step
time.steps.m_delta = 10000;              % Size of big time steps: 240fs
time.steps.s_number = 400;               % Number of sub steps per time step

% Hamiltonian operator 
hamilt.truncate.e_min = -2.0;            % Lower truncation of energy
hamilt.truncate.e_max = +1.5;            % Upper truncation of energy

hamilt.pot{1,1}    = pot.metiu1;         % Potentials: Cylindrical coordinates
hamilt.pot{1,1}.a0 = +700/atomic.w.cm_1; % constant offset
hamilt.pot{1,1}.a1 = +100/atomic.w.cm_1; % first coefficient
hamilt.pot{1,1}.a2 = -600/atomic.w.cm_1; % second coefficient

% For qm_bound only
hamilt.eigen.stop = 15;

% Initial wave function is a spherical harmonic of degree l=1.
time.dof{1}        = init.fbr;           % use eigenstate of grid (spherical harmonic)
time.dof{1}.state  = 2;                  % use second eigenstate (for m=0, this is l=1)
