% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global atomic hamilt plots space time

prt.disp ( '*************************************' )
prt.disp ( 'HCl+ cation ( X==>A excitation )' )
prt.disp ( 'A.D.Pradhan, K.P.Kirby, A.Dalgarno  )' )
prt.disp ( 'J. Chem. Phys. 95, 9009 (1991)      )' )
prt.disp ( '*************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs      = 2;
hamilt.coupling.represent  = 'dia';
hamilt.coupling.labels     = {'X ^2\Pi', 'A ^2\Sigma^+'};
hamilt.coupling.ini_coeffs = [0 1 0 0 0 0 0];   % Initially: ground state only
hamilt.coupling.ini_rep    = 'dia';

% Spatial discretization
space.dof{1} = dof.fft;                  % using FFT grid
space.dof{1}.mass = 0.97989/atomic.m.u;  % Reduced mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min =  1.5;               % Lower bound of grid 
space.dof{1}.x_max =  7.5;               % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 040;               % Index of final time step
time.steps.m_delta  = 25/atomic.t.as;    % Size of time steps: 25 as 
time.steps.s_number =   025;             % Number of sub steps per time step

% Electric field as sequence of pulses
time.efield.dressed = true;              % Using dressed state picture
time.efield.photons = {-2:2:2, -3:2:3};  % Number of photons for 3+4=7 dressed states
time.pulse{1}       = efi.sin_2;         % Shape of envelope
time.pulse{1}.delay = 0.50/atomic.t.fs;  % Time delay of pulse center: 0.25 fs
time.pulse{1}.fwhm  = 0.50/atomic.t.fs;  % Full width at half maximum: 0.25 fs
time.pulse{1}.ampli = 0.5;               % Field amplitude
time.pulse{1}.frequ = 3.52/atomic.E.eV;  % Carrier frequency: 3.52 eV
time.pulse{1}.phase = +time.pulse{1}.delay*time.pulse{1}.frequ+pi/2;    % Phase

% Hamiltonian operator 
hamilt.truncate.e_min    =  -2.5;        % Lower truncation of energy
hamilt.truncate.e_max    =  -1.0;        % Upper truncation of energy

hamilt.pot{1,1}        = pot.interp;     % interpolate tabulated potential
hamilt.pot{2,2}        = pot.interp;     % interpolate tabulated potential

hamilt.dip{1}          = cell(2);
hamilt.dip{1}{1,2}     = dip.interp;     % Use tabulated values

% Initial wave function
time.dof{1}        = init.gauss;         % Gaussian-shaped wavepacket
time.dof{1}.width  =  0.153688;          % Width 
time.dof{1}.pos_0  =  2.519868;          % Center in position representation
time.dof{1}.mom_0  =  0.0;               % Center in momentum representation

% Plot densities
plots.density        = vis.contour;
plots.density.expect = false;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_min = -2;
plots.expect.e_max = 0.5;
