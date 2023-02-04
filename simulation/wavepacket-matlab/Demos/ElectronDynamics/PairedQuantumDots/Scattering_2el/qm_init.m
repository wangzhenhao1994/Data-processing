% Copyright (C) 2020 Federico Pont
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init()
global hamilt plots space time
 
prt.disp ( '***************************************************************' )
prt.disp ( 'Two electrons in a Paired Quantum Dot with e-e interaction' )
prt.disp ( '***************************************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % Fourier grid
space.dof{1}.n_pts = 384;                % Number of grid points
space.dof{1}.x_min = -300;               % Lower bound of the grid
space.dof{1}.x_max =  300;               % Upper bound of the grid
space.dof{1}.mass  =  1;                 % Mass for the kinetic energy

space.dof{2}       = dof.fft;            % Fourier grid
space.dof{2}.n_pts = 384;                % Number of grid points
space.dof{2}.x_min = -300;               % Lower bound of the grid
space.dof{2}.x_max =  300;               % Upper bound of the grid
space.dof{2}.mass  = 1;                  % Mass for the kinetic energy

% Hamiltonian operator: potential energy function
hamilt.pot{1,1} = pot.pqd1D;             % Paired Quantum dots expansion
hamilt.pot{1,1}.b_L = 0.3;               % Inverse square Size of the left dot
hamilt.pot{1,1}.b_R = 1.0;               % Inverse square Size of the right dot
hamilt.pot{1,1}.V_L = -0.65;             % Depth of right quantum dot
hamilt.pot{1,1}.V_R = -0.6;              % Depth of left quantum dot
hamilt.pot{1,1}.R = 10;                  % Distance between the quantum dot centers
hamilt.pot{1,1}.lambda = 1;              % Coulomb interaction

% Absorbing boundary conditions
hamilt.nip{1}  = nip.power;              
hamilt.nip{1}.exp  = [2 2];
hamilt.nip{1}.min  = [-200 -200];
hamilt.nip{1}.max  = [200 200];

% Additional multiplicative operator: Reaction probabilities
hamilt.amo{1} = amo.reaction;            % population of reactant states
hamilt.amo{1}.reac = 2;                  % educt distance is DOF #2
hamilt.amo{1}.prod = 1;                  % product distance is DOf #1
hamilt.amo{1}.side = 'r';                % projection on the reactants
hamilt.amo{1}.label= 'reflected';        % labeling the reactants

hamilt.amo{2} = amo.reaction;            % population of product states
hamilt.amo{2}.reac = 2;                  % educt distance is DOF #2
hamilt.amo{2}.prod = 1;                  % product distance is DOf #1
hamilt.amo{2}.side = 'p';                % projection on the products
hamilt.amo{2}.label= 'transmitted';      % labeling the products

hamilt.amo{3} = amo.interval;            % population of product states
hamilt.amo{3}.min = [-10 10];            % educt distance is DOF #2
hamilt.amo{3}.max = [0 200];             % product distance is DOf #1
hamilt.amo{3}.label= 'ICEC1';            % labeling the products

hamilt.amo{4} = amo.interval;            % population of product states
hamilt.amo{4}.min = [10 -10];            % educt distance is DOF #2
hamilt.amo{4}.max = [200 0];             % product distance is DOf #1
hamilt.amo{4}.label= 'ICEC2';            % labeling the products

% Temporal discretization
time.steps.m_start    = 000;             % first time step
time.steps.m_stop     = 030;             % last time step
time.steps.m_delta    = 025;             % length of one time step
time.steps.s_number   = 500;             % number of sub steps

% Initial wave function: Superposition of two Gaussians
time.corr        = init.gauss_corr;
time.corr.pos_0  = [-125,5; 5, -125];
time.corr.mom_0  = [.51,0; 0, .51];
time.corr.width  = [10, 1.18; 1.18, 10];
time.corr.factor = [1 +1];               % symmetric or anti-symmetric

% Plots of time-evolving densities
plots.density          = vis.contour;    % contour or surface plots
plots.density.represent = 'dvr';         % Position (dvr) or momentum (fbr) densities
plots.density.complex = 'abs2';          % Choice of real|imag|abs2
plots.density.cnt_nlev = [60 15];
plots.density.pot_max = 2;

% Plot of expectation values
plots.expect   = vis.expect;
plots.expect.e_min = -1;
plots.expect.e_max = +1;
plots.expect.p_min = -0.001;
plots.expect.p_max = +0.010;
