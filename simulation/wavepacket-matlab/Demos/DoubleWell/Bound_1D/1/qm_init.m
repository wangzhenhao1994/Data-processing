% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it 
% into other works without any restrictions.

function qm_init
global hamilt plots space

prt.disp ( '**************************************' )
prt.disp ( 'Razavy symmetric double well potential' )
prt.disp ( 'Reference calculation using FFT-DVR' )
prt.disp ( '**************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % using FFT grid
space.dof{1}.mass  = 1/2;                % Particle mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -7.0;               % Lower bound of grid 
space.dof{1}.x_max =  7.0;               % Upper bound of grid

% Hamiltonian operator 
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
% hamilt.eigen.symmetry = 'g';             % Even parity only

% Plot of densities
plots.density         = vis.contour;
plots.density.pot_max = 100;
plots.density.expect  = false;

% Plot of expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 80;
