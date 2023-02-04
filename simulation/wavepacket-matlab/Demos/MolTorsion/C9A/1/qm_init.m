% Copyright (C) 2009,2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init(symmetry)
global atomic hamilt plots space

prt.disp ( '**********************************************' )
prt.disp ( 'Calculation of the vibrational spectrum of the' )
prt.disp ( 'S0 state of 9-(N-carbazolyl)-anthracene (C9A) ' )
prt.disp ( 'see Z.Phys.D 34:111                           ' )
prt.disp ( '**********************************************' )

% moment of intertia is hbar/(4pi*c * 0.035 cm^-1)
c = 299792458/atomic.l.m*atomic.t.s; % approx 137
inertia = 1/(4 * pi * c * 0.035 * atomic.l.A * 1e-8);

% Spatial discretization
space.dof{1}       = dof.fft;            % Fourier grid
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = 030*pi/180;         % Lower bound of grid 30 degrees
space.dof{1}.x_max = 150*pi/180;         % Upper bound of grid 150 degrees
space.dof{1}.mass  = inertia;            % Moment of inertia for the kinetic energy

% Hamiltonian operator 
hamilt.truncate.e_min    = -1.0;         % Lower truncation of energy
hamilt.truncate.e_max    = 1.0;          % Upper truncation of energy

hamilt.pot{1,1} = pot.C9A ('S0');        % ground state surface of C9A
hamilt.pot{1,1}.V0 =    +17.00 / atomic.w.cm_1; 
hamilt.pot{1,1}.V2 =  -1430.86 / atomic.w.cm_1;
hamilt.pot{1,1}.V4 = 180648.43 / atomic.w.cm_1;

% Select eigen/values/functions
hamilt.eigen.start     =  0;             % Lower index
hamilt.eigen.stop      =  0;             % Upper index
hamilt.eigen.symmetry  = symmetry;       % get symmetry passed in by global variable

% Plot densities
plots.density        = vis.contour;     
plots.density.expect = false;

% Plot expectation values
plots.expect = vis.expect;