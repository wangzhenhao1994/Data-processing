% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots space time

prt.disp ( '*********************************************' )
prt.disp ( 'Tabulated potential (neutral state of Cl-NH3)' )
prt.disp ( 'Physical Review Letters 93 (4), 048301 (2004)' )
prt.disp ( '*********************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % Fourier grid
space.dof{1}.label = 'NH_2-H';           % Axis label
space.dof{1}.x_min = 0.25/atomic.l.A;    % Lower bound of grid
space.dof{1}.x_max = 3.00/atomic.l.A;    % Upper boundof grid
space.dof{1}.n_pts = 128;                % Number of grid points

space.dof{2}       = dof.fft;            % Fourier grid
space.dof{2}.label = 'H-Cl';             % Axis label
space.dof{2}.x_min = 0.70/atomic.l.A;    % Lower bound of grid
space.dof{2}.x_max = 3.10/atomic.l.A;    % Upper boundof grid
space.dof{2}.n_pts = 128;                % Number of grid points

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 020;               % Index of final time step
time.steps.m_delta  = 050;               % Size of time steps 
time.steps.s_number = 100;               % Number of sub steps per time step

% Hamiltonian operator 
hamilt.truncate.e_min    =  -0.05;       % Lower truncation of energy
hamilt.truncate.e_max    =  +0.20;       % Upper truncation of energy

hamilt.kin{1}       = kin.triatomic;     % Jacobi coordinates with fixed angle
hamilt.kin{1}.dof   = [1 2];             % on which degrees of freedom we act
hamilt.kin{1}.mass  = [16 1 35] / atomic.m.u; % masses of NH2,H,Cl
hamilt.kin{1}.theta = pi;                % NH2-H-Cl bending angle

hamilt.pot{1,1}          = pot.interp;   % Interpolation of potential fct.
hamilt.pot{1,1}.n_pts    = [56 49];      % Number of tabulated values
hamilt.pot{1,1}.pos_conv = atomic.l.A;   % Conversion of coordinates
hamilt.pot{1,1}.pot_conv = atomic.w.cm_1;% Conversion of pot. energies

% Initial wave function
time.dof{1}        = init.gauss;
time.dof{1}.width  = 0.28 / 2;
time.dof{1}.pos_0  = 1.95;
time.dof{1}.mom_0  = 0;

time.dof{2}        = init.gauss;
time.dof{2}.width  = 0.47 / 2;
time.dof{2}.pos_0  = 4.52;
time.dof{2}.mom_0  = 0;

% Plot densities
plots.density          = vis.contour;    % Choose plot type
% plots.density.srf_view = [60 65];        % View point for surface plot: [az el]

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 0.015;              % Manually set energy range 
