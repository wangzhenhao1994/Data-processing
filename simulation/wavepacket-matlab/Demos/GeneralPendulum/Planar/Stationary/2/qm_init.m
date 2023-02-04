% Copyright (C) 2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space

prt.disp ( '***************************************' )
prt.disp ( 'Supersymmetry and eigensurface topology' )
prt.disp ( 'of the planar quantum pendulum' )
prt.disp ( 'B. Schmidt and B. Friedrich' )
prt.disp ( 'Front. Phys. 2, 37' )
prt.disp ( 'DOI:10.3389/fphy.2014.00037' )
prt.disp ( 'Reproducing red circles in Fig. 4' )
prt.disp ( 'USING HERE: FFT-DVR for THETA' )
prt.disp ( '***************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % Fourier grid
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -2*pi;              % Lower bound of grid
space.dof{1}.x_max = +2*pi;              % Upper bound of grid
space.dof{1}.mass  =  1/2;               % Mass for the kinetic energy

% Orientation
hamilt.amo{1} = amo.cosine;              % cosine projector
hamilt.amo{1}.exp = 1;                   % exponent

% Alignment
hamilt.amo{2} = amo.cosine;              % cosine^2 projector
hamilt.amo{2}.exp = 2;                   % exponent

% Potential: beta=5, kappa=6
hamilt.pot{1,1} = pot.pendulum;          % Potential for generalized pendula
hamilt.pot{1,1}.eta  = 30;               % Orientation: cos
hamilt.pot{1,1}.zeta = 25;               % Alignment: cos^2

% Generate first 25 eigenstates
hamilt.eigen.start = 0;
hamilt.eigen.stop = 25;

% Plot time evolution of density
plots.density         = vis.contour;     % Wigner contour plot
plots.density.pot_max = 050;             % customize density plot
plots.density.expect  = false;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 040;                % Set maximum for energy plot
