% Copyright (C) 2009,2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init(symmetry)
global atomic hamilt plots space

prt.disp ( '********************************************************************' )
prt.disp ( 'Calculation of absorption spectrum of the S1 state only' )
prt.disp ( 'of 9-(N-carbazolyl)-anthracene (C9A)                   ' )
prt.disp ( 'see Z.Phys.D 34:111                                    ' )
prt.disp ( '********************************************************************' )

% moment of intertia is hbar/(4pi*c * 0.035 cm^-1)
c = 299792458/atomic.l.m*atomic.t.s; % approx 137
inertia = 1/(4 * pi * c * 0.035 * atomic.l.A * 1e-8);

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs = 2;
hamilt.coupling.represent = 'dia';

% Spatial discretization
space.dof{1}       = dof.fft;            % Fourier grid
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = 030*pi/180;         % Lower bound of grid 30 degrees
space.dof{1}.x_max = 150*pi/180;         % Upper bound of grid 150 degrees
space.dof{1}.mass  = inertia;            % Moment of inertia for the kinetic energy

% Hamiltonian operator 
hamilt.truncate.e_min    = -1.0;         % Lower truncation of energy
hamilt.truncate.e_max    = +1.0;         % Upper truncation of energy

hamilt.pot{1,1}        = pot.C9A ('S1'); % First excited state of C9A
hamilt.pot{1,1}.coeffs = [ 6712 -2777 -6963 -2007 1692 373 ...
    -4688 -7518 -5501 -602 3348 4539 3511 1873 662 140] / atomic.w.cm_1;
hamilt.pot{1,1}.vshift = (- 193.3 + 25893.3) / atomic.w.cm_1;

hamilt.pot{2,2}      = pot.C9A ('Sx');   % "dark" excited state of C9A
hamilt.pot{2,2}.V0   = 25806 / atomic.w.cm_1;
hamilt.pot{2,2}.V2   = 19094 / atomic.w.cm_1;

hamilt.pot{1,2}      = pot.C9A ('S1-Sx');% S1-Sx coupling of C9A
hamilt.pot{1,2}.C    = 22/atomic.w.cm_1; % coupling strength in atomic units
hamilt.pot{1,2}.a    = 0.056;            % width of Gaussian
hamilt.pot{1,2}.phi0 = 0.227;            % center of Gaussian
                     
% Select eigen/values/functions
hamilt.eigen.start     =  0;             % Lower index
hamilt.eigen.stop      =  35;            % Upper index
hamilt.eigen.symmetry  = symmetry;       % Take the symmetry from a global variable.

% Plot densities
plots.density        = vis.contour;     
plots.density.expect = false;

% Plot expectation values
plots.expect = vis.expect;
plots.expect.e_min = 0.115;
plots.expect.e_max = 0.125;
