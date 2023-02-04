% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space time

prt.disp ( '*********************************************************************' )
prt.disp ( 'Evolution of pendulum with squeezed initial function    ' )
prt.disp ( '*********************************************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % using fft grid
space.dof{1}.mass  = 1;                  % (Reduced) moment of inertia
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = 0;                  % Lower bound of grid 
space.dof{1}.x_max = 2*pi;               % Upper bound of grid

% Additional multiplicative operator
hamilt.amo{1}     = amo.cosine;          % Use cos-function as AMO
hamilt.amo{1}.exp = 1;                   % Degree of alignment

% Temporal discretization
time.steps.m_start = 000;                % Index of initial time step
time.steps.m_stop  = 100;                % Index of final time step
time.steps.m_delta = pi/50;              % Size of time steps
time.steps.s_number = 05;                % Number of sub steps per time step

% Hamiltonian operator 
hamilt.truncate.e_min  = 000;            % Lower truncation of energy
hamilt.truncate.e_max  = 100;            % Upper truncation of energy

hamilt.pot{1,1}     = pot.pendulum;      % Plane pendulum
hamilt.pot{1,1}.eta = -5.0;              % Prefactor of cos-potential
hamilt.pot{1,1}.v_0 = +5.0;              % Energy offset

% Initial wave function
time.dof{1}          = init.pendulum1;   % Mathieu type wavefunction
time.dof{1}.parity   = 'c';              % cosine elliptic  
time.dof{1}.order    = 0;                % Order of Mathieu function 
time.dof{1}.multiple = 1;                % Potential multiplicity
time.dof{1}.barrier  = 100;              % Potential barrier
time.dof{1}.shift    = 0;                % Potential shift (horizontal)

% Plot densities
plots.density      = vis.polar;          % Polar plot 

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 003;                % Range for energy plot

