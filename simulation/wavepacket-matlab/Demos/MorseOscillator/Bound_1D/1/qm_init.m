% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space

prt.disp ( '*****************************' )
prt.disp ( 'Morse oscillator (OH radical)' )
prt.disp ( '*****************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % using fft grid
space.dof{1}.mass = 1728.539;            % Reduced mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min =  1.0;               % Lower bound of grid 
space.dof{1}.x_max = 10.0;               % Upper bound of grid
space.dof{1}.periodic = false;           % Build the kinetic energy matrix
                                         % without periodic boundary conditions
% Hamiltonian operator 
hamilt.truncate.e_min  =  0.0;           % Lower truncation of energy
hamilt.truncate.e_max  =  1.0;           % Upper truncation of energy

hamilt.pot{1,1}      = pot.morse;        % Harmonic oscillator
hamilt.pot{1,1}.d_e  = 0.1994;           % Dissociation energy
hamilt.pot{1,1}.r_e  = 1.821;            % Equilibrium length
hamilt.pot{1,1}.alf  = 1.189;            % Range parameter

% Select eigen/values/functions
hamilt.eigen.start     = 00;             % Lower index
hamilt.eigen.stop      = 21;             % Upper index

% Plot time evolution of the density
plots.density         = vis.contour;     % Contour plot of Wigner function
plots.density.expect  = false;

% Plot time evolution of expectation values
plots.expect          = vis.expect;
plots.expect.e_max    = 0.25;            % Manually set range for energy plot
