% Copyright (C) 2020 Federico Pont
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init()
global hamilt plots space 

prt.disp ( '***************************************************************' )
prt.disp ( 'One electron in a Paired Quantum Dot ' )
prt.disp ( '***************************************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % Fourier grid
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = -20;                % Lower bound of the grid
space.dof{1}.x_max =  20;                % Upper bound of the grid
space.dof{1}.mass  =  1;                 % Mass for the kinetic energy

% Hamiltonian operator: potential energy function
hamilt.pot{1,1} = pot.pqd1D;             % Paired Quantum dots potential
hamilt.pot{1,1}.b_L = 0.3;               % Inverse square Size of the left dot
hamilt.pot{1,1}.b_R = 1.0;               % Inverse square Size of the right dot
hamilt.pot{1,1}.V_L = -0.71;             % Depth of left quantum dot
hamilt.pot{1,1}.V_R = -0.60;             % Depth of right quantum dot
hamilt.pot{1,1}.R = 10;                  % Distance between the quantum dot centers

% Select eigen/values/functions
hamilt.eigen.start    = 0;               % Lower index
hamilt.eigen.stop     = 20;              % Upper index

% Plots of bound state densities
plots.density          = vis.curve;    % 3D surface plot
plots.density.pot_max = 1.5;

% Plot of expectation values
plots.expect   = vis.expect;
plots.expect.e_max = 2;
plots.expect.e_min = -1;
