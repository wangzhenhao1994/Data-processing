% Copyright (C) 2020 Federico Pont
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init()
global hamilt plots space

prt.disp ( '***************************************************************' )
prt.disp ( 'Two electrons in a Paired Quantum Dot with e-e interaction' )
prt.disp ( '***************************************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % Fourier grid
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = -20;                % Lower bound of the grid
space.dof{1}.x_max =  20;                % Upper bound of the grid
space.dof{1}.mass  =  1;                 % Mass for the kinetic energy

space.dof{2}       = dof.fft;            % Fourier grid
space.dof{2}.n_pts = space.dof{1}.n_pts; % Number of grid points
space.dof{2}.x_min = space.dof{1}.x_min; % Lower bound of the grid
space.dof{2}.x_max = space.dof{1}.x_max; % Upper bound of the grid
space.dof{2}.mass  = space.dof{1}.mass;  % Mass for the kinetic energy

% Hamiltonian operator: potential energy function
hamilt.pot{1,1} = pot.pqd1D;             % Paired Quantum dots expansion
hamilt.pot{1,1}.b_L = 0.3;               % Inverse square Size of the left dot
hamilt.pot{1,1}.b_R = 1.0;               % Inverse square Size of the right dot
hamilt.pot{1,1}.V_L = -0.71;             % Depth of right quantum dot
hamilt.pot{1,1}.V_R = -0.6;              % Depth of left quantum dot
hamilt.pot{1,1}.R = 10;                  % Distance between the quantum dot centers
hamilt.pot{1,1}.lambda = 1;              % Coulomb interaction

% Select eigen/values/functions
hamilt.eigen.start    = 000;             % Lower index
hamilt.eigen.stop     = 020;             % Upper index

% Plots of bound state densities
plots.density          = vis.contour;
plots.density.pot_max = 2;

% Plot of expectation values
plots.expect   = vis.expect;
plots.expect.e_max = 0;
plots.expect.e_min = -1;
