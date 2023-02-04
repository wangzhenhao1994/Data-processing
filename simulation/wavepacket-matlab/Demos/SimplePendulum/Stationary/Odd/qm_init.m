% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space

prt.disp ( '***************************************' )
prt.disp ( 'Eigenstates of plane quantum pendulum  ' )
prt.disp ( '***************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % using fft grid
space.dof{1}.mass  = 1;                  % mass
space.dof{1}.n_pts = 192;                % Number of grid points
space.dof{1}.x_min = 0;                  % Lower bound of grid 
space.dof{1}.x_max = 2*pi;               % Upper bound of grid

hamilt.amo{1}     = amo.cosine;          % Use cos-function as AMO
hamilt.amo{1}.exp = 1;                   % Degree of alignment

% Hamiltonian operator 
hamilt.truncate.e_min = 000;             % Lower truncation of energy
hamilt.truncate.e_max = 300;             % Upper truncation of energy

hamilt.pot{1,1}     = pot.pendulum;      % Intramolecular torsion
hamilt.pot{1,1}.eta = -50;               % Prefactor of cos-potential
hamilt.pot{1,1}.v_0 = +50;               % Energy offset

% Variables concerning the eigenvalue problem: Select eigen/values/functions
hamilt.eigen.start     =  0;             % Lower index
hamilt.eigen.stop      = 20;             % Upper index
hamilt.eigen.symmetry  = 'u';            % Odd parity states

% Plot densities
plots.density          = vis.polar;      % Create object

% Plot expectation values
plots.expect           = vis.expect;     % Create object
plots.expect.e_max     = 300;
