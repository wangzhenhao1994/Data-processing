% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space

prt.disp ( '**********************************************' )
prt.disp ( 'Razavy symmetric double well potential' )
prt.disp ( 'Reference calculation using Gauss-Hermite DVR ' )
prt.disp ( '**********************************************' )

% Spatial discretization
space.dof{1}       = dof.hermite;        % Use Gauss-Hermite grid
space.dof{1}.mass  = 1/2;                % Particle mass
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.v_2   = 7;                  % Frequency of oscillator

% Hamiltonian operator 
% hamilt.eigen.symmetry = 'g';             % Symmetry may be exploited
hamilt.truncate.e_min  = -15.0;          % Lower truncation of energy
hamilt.truncate.e_max  = 1000.0;         % Upper truncation of energy

% Razavy Potential: beta=0.1, kappa=-7
hamilt.pot{1,1}          = pot.razavy;   % Hyperbolic potential
hamilt.pot{1,1}.modified = true;         % Use modified version
hamilt.pot{1,1}.eta      = -0.7;         % prefactor of cosh
hamilt.pot{1,1}.zeta     = 0.01;         % prefactor of cosh^2

% Select eigen/values/functions
hamilt.eigen.start    =  0;              % Lower index
hamilt.eigen.stop     = 30;              % Upper index

% Plot time evolution of density
plots.density         = vis.curve;       % Colored curve plot
plots.density.pot_max = 100;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 80;
