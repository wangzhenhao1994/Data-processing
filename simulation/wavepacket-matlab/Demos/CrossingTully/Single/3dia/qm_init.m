% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~)
global hamilt plots space time

prt.disp ( '****************************************************************' )
prt.disp ( 'Single crossing example, k=10 (Tully 1990)' )
prt.disp ( '****************************************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs      = 2;          % Two coupled channels
hamilt.coupling.represent  = 'dia';      % Choose diabatic representation
hamilt.coupling.ini_rep    = 'dia';      % Choose diabatic initial state
hamilt.coupling.ini_coeffs = [1 0];      % Starting on lower left diabatic state 

% Spatial discretization
space.dof{1} = dof.fft;                  % Using FFT based DVR/FBR scheme
space.dof{1}.mass = 2000;                % Mass of particle
space.dof{1}.n_pts = 384;                % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max =  20.0;              % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 060;               % Index of final time step
time.steps.m_delta  = 060;               % Size of time steps 
time.steps.s_number = 030;               % Number of sub steps per time step

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width =  0.50;               % Width 
time.dof{1}.pos_0 = -6.0;                % Center in position representation
time.dof{1}.mom_0 = 10.0;                % Center in momentum representation

% Hamiltonian operator 
hamilt.truncate.e_min    = -0.01;        % Lower truncation of energy
hamilt.truncate.e_max    =  0.10;        % Upper truncation of energy

hamilt.pot{1,1} = pot.tully1;            % Single Crossing
hamilt.pot{2,2} = pot.tully1;            % Single Crossing
hamilt.pot{1,2} = pot.tully1;            % Single Crossing

% Absorbing boundary conditions
for m=1:2                                % Same for both channels
    hamilt.nip{m}      = nip.power;      % Negative imaginary potential
    hamilt.nip{m}.exp  = 3;              % Exponent
    hamilt.nip{m}.min = -08.0;           % Beginning of inner grid region
    hamilt.nip{m}.max = +18.0;           % End of inner grid region
end

% Plot densities
plots.density           = vis.curve;     % Contour plot of the Wigner function
plots.density.represent = 'dvr';         % Position or momentum space
plots.density.complex   = 'abs2';        % Real or Imag or Abs2
plots.density.energy    = true;          % Show also energy function
plots.density.range     = false;         % Manual setting of plotting range
plots.density.pot_max   = 0.050;

% plots.density         = vis.contour;
% plots.density.expect  = false;
% plots.density.range   = true;             % manual setting of plotting range
% plots.density.x_min   = -10;
% plots.density.x_max   = 10;
% plots.density.y_min   = 0;
% plots.density.y_max   = 15;
% plots.density.pot_max = 0.015;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 0.04;
