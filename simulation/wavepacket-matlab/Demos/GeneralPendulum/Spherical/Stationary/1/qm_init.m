% Copyright (C) 2016 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (b,m)
global hamilt plots space

prt.disp ( '***************************************' )
prt.disp ( 'Supersymmetry and eigensurface topology' )
prt.disp ( 'of the spherical quantum pendulum' )
prt.disp ( 'B. Schmidt and B. Friedrich' )
prt.disp ( 'Phys. Rev. A 91, 022111' )
prt.disp ( 'DOI:10.1103/PhysRev.A.91.022111' )
prt.disp ( 'Reproducing red circles in Fig. 3' )
prt.disp ( 'USING HERE: FFT-DVR for THETA' )
prt.disp ( '***************************************' )

% Spatial discretization
space.dof{1}       = dof.fft;            % Fourier grid
space.dof{1}.label = '\Theta';
space.dof{1}.n_pts = 1024;               % Number of grid points
space.dof{1}.x_min = 0;                  % Lower bound of grid
space.dof{1}.x_max = +1*pi;              % Upper bound of grid
space.dof{1}.mass  =  1/2;               % Mass for the kinetic energy

% Orientation
hamilt.amo{1} = amo.cosine;              % cosine projector
hamilt.amo{1}.exp = 1;                   % exponent

% Alignment
hamilt.amo{2} = amo.cosine;              % cosine^2 projector
hamilt.amo{2}.exp = 2;                   % exponent

% Pendular potential
hamilt.pot{1,1} = pot.pendulum;          % Potential for generalized pendula
hamilt.pot{1,1}.xi = m^2 - 1/4;          % Azimuthal rotation
hamilt.pot{1,1}.eta = 2*b*(m+1);         % Orientation: cos Theta
hamilt.pot{1,1}.zeta = b^2;              % Alignment: cos^2 Theta
hamilt.pot{1,1}.v_0 = -1/4;              % Energy shift

% Plot time evolution of density
plots.density         = vis.curve;       % Colored curve plot
plots.density.pot_max = 200;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.p_min = -0.3;               % customize population/amo plot
plots.expect.p_max = +1.1;               % customize population/amo plot
plots.expect.e_min = -150;               % customize energy plot
plots.expect.e_max = +150;               % customize energy plot


