% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space

prt.disp ( '************************************' )
prt.disp ( 'Harmonic oscillator in one dimension' )
prt.disp ( '************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % Fourier grid
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = -7;                 % Lower bound of grid
space.dof{1}.x_max =  7;                 % Upper bound of grid
space.dof{1}.mass  =  1;                 % Mass for the kinetic energy

hamilt.truncate.e_min  =  0.0;           % Lower truncation of energy
hamilt.truncate.e_max  = 50.0;           % Upper truncation of energy

hamilt.pot{1,1}        = pot.taylor;     % Taylor expansion
hamilt.pot{1,1}.coeffs = [0;1];          % Force constant

% Select eigen/values/functions, optional symmetry
hamilt.eigen.start        =  0;          % Lower index
hamilt.eigen.stop         = 10;          % Upper index
hamilt.eigen.symmetry     = 'n';         % Enforce symmetry

% Modify settings for appearance of plots (if desired)
plots.density.surface.view  = [30 60];   % View point for surface plot: [az el]

% Plots of densities
plots.density         = vis.contour;     % contour plot
plots.density.pot_max = 15;              % Tweaking of density display
plots.density.expect  = false;           % Do not draw expectation values

% Plots of expectation values
plots.expect        = vis.expect;   
plots.expect.e_max  = 11;                % Tune the display of the energy values
