% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space

prt.disp ( '*************************************' )
prt.disp ( 'Harmonic oscillator in two dimensions' )
prt.disp ( '*************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % Fourier grid
space.dof{1}.n_pts = 48;                 % Number of grid points
space.dof{1}.x_min = -7;                 % Lower bound of the grid
space.dof{1}.x_max =  7;                 % Upper bound of the grid
space.dof{1}.mass  =  1;                 % Mass for the kinetic energy

space.dof{2}       = dof.fft;            % Fourier grid
space.dof{2}.n_pts = 48;                 % Number of grid points
space.dof{2}.x_min = -7;                 % Lower bound of the grid
space.dof{2}.x_max =  7;                 % Upper bound of the grid
space.dof{2}.mass  =  1;                 % Mass for the kinetic energy

% Hamiltonian operator 
hamilt.truncate.e_min  =  0.0;           % Lower truncation of energy
hamilt.truncate.e_max  = 50.0;           % Upper truncation of energy

hamilt.pot{1,1}        = pot.taylor;     % Taylor expansion
hamilt.pot{1,1}.coeffs = [0 0; 0.81 1];  % Force constants

% Select eigen/values/functions, optional symmetry
hamilt.eigen.start    = 000;             % Lower index
hamilt.eigen.stop     = 025;             % Upper index
hamilt.eigen.symmetry = 'n';             % Even parity only
hamilt.eigen.storage   = 's';            % sparse storage

% Plots of densities
plots.density          = vis.surface;    % 3D surface plot
plots.density.srf_view = [25 55];        % View point for surface plot: [az el]

% Plots of expectation values
plots.expect        = vis.expect;
plots.expect.e_max  = 08;                % manually set range for energy plot
