% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008 Burkhard Schmidt
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots space

prt.disp ( '*********************************************' )
prt.disp ( 'Fulvene torsion: Stationary torsion (A state)' )
prt.disp ( 'See http://dx.doi.org/10.1002/cphc.200600543 ' )
prt.disp ( '*********************************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % using fft grid
space.dof{1}.mass = 1/(2.44/atomic.E.meV);  % Reduced moment of inertia: 1/2.44 meV
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -pi;                % Lower bound of grid 
space.dof{1}.x_max =  pi;                % Upper bound of grid

% Hamiltonian operator 
hamilt.truncate.e_min    =  0.0;         % Lower truncation of energy
hamilt.truncate.e_max    = 0.25;         % Upper truncation of energy

hamilt.pot{1,1}    = pot.pendulum;       % Intramolecular torsion
hamilt.pot{1,1}.zeta = - 856*2.44/atomic.E.meV; % Barrier height: 856*2.44 meV

% Select eigen/values/functions
hamilt.eigen.start     = 0;              % Lower index
hamilt.eigen.stop      =10;              % Upper index
hamilt.eigen.symmetry  = 'u';            % Odd parity eigenstates

% Plot densities
plots.density        = vis.contour;      % Contour plot of the eigenstates
plots.density.expect = false;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 0.05;               % manually set the ranges for the energy plot

