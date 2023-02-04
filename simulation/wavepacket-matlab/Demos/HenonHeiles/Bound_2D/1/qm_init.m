% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space

prt.disp ( '*************************************' )
prt.disp ( 'Henon-Heiles system in two dimensions' )
prt.disp ( 'see, e.g., M.J.Davis and E.J.Heller  ' )
prt.disp ( 'J. Chem. Phys. 71, 3383-3395 (1979)  ' )
prt.disp ( '*************************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % using fft grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 48;                 % Number of grid points
space.dof{1}.x_min = -8.0;               % Lower bound of grid 
space.dof{1}.x_max = +8.0;               % Upper bound of grid
space.dof{1}.periodic = false;           % use non-periodic kinetic energy

% Spatial discretization
space.dof{2} = dof.fft;                  % using fft grid
space.dof{2}.mass = 1;                   % mass
space.dof{2}.n_pts = 48;                 % Number of grid points
space.dof{2}.x_min = -6.0;               % Lower bound of grid 
space.dof{2}.x_max = +9.0;               % Upper bound of grid
space.dof{1}.periodic = false;           % use non-periodic kinetic energy

% Hamiltonian operator 
hamilt.truncate.e_min  =   0.0;          % Lower truncation of energy
hamilt.truncate.e_max  = +30.0;          % Upper truncation of energy

hamilt.pot{1,1}   = pot.henon;           % Henon-Heiles system
hamilt.pot{1,1}.A = 1;                   % Force constant
hamilt.pot{1,1}.L = 1/(4*sqrt(5));       % Cubic parameter

% Select eigen/values/functions
hamilt.eigen.start    = 000;             % Lower index
hamilt.eigen.stop     = 098;             % Upper index

% Plotting densities
plots.density           = vis.surface;   % draw it as a 3D plot
plots.density.srf_view  = [60 55];       % View point for surface plot: [az el]

% Plotting expectation values 
plots.expect        = vis.expect;
plots.expect.e_max  = 15;