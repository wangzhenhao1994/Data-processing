% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%                    Boris Schaefer-Bung
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (suffix)
global control hamilt plots space state time reduce

if nargin>0
    state.save_suffix = suffix;
end

prt.disp ( '**********************************' )
prt.disp ( 'Asymmetric double well potential  ' )
prt.disp ( ' see example 3 in:                ' )
prt.disp ( ' B. Schaefer-Bung, C. Hartmann,   ' )
prt.disp ( ' B. Schmidt, and Ch. Schuette,    ' )
prt.disp ( ' J. Chem. Phys. 135, 014112 (2011)' )
prt.disp ( '**********************************' )

% Spatial discretization
space.dof{1} = dof.fft;                  % using FFT grid
space.dof{1}.mass = 162;                 % I=162 [hbar^2/D], see section IV.D
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -2.0;               % Lower bound of grid 
space.dof{1}.x_max =  2.0;               % Upper bound of grid

% Temporal discretization
time.steps.m_start = 0;                  % Index of initial time step
time.steps.m_stop  = 100;                % Index of final time step
time.steps.m_delta = 7.2;                % Size of time steps 

% Electric field as half-cycle pulse
time.pulse{1}       = efi.sin_2;         % sin^2 shape of envelope
time.pulse{1}.delay = 180;               % Time delay of pulse center
time.pulse{1}.fwhm  = 180;               % Full width at half maximum
time.pulse{1}.ampli = 0.2;               % u_0=0.2 [D/(e*a_0)] see section IV.D
time.pulse{1}.frequ = 0.196;             % Omega = 0.196 [D/hbar], see section IV.D
time.pulse{1}.phase = 35.28;             % Phase = frequ * delay = 0.196* 180

% Hamiltonian operator 
hamilt.truncate.e_min  =  -1.0;          % Lower truncation of energy
hamilt.truncate.e_max  =   10.0;         % Upper truncation of energy

hamilt.pot{1,1}        = pot.taylor;     % Taylor series: Double well potential
hamilt.pot{1,1}.coeffs = [0.055;-4;0;+24];  % Linear, quadratic, quartic constant
hamilt.pot{1,1}.vshift = 1;              % Constant energy shift

hamilt.dip{1}{1,1}        = dip.taylor;  % Dipole moment: Taylor series
hamilt.dip{1}{1,1}.coeffs = 1;           % Linear dipole moment, slope 1

hamilt.sbc{1,1}        = sbc.taylor;     % System-bath coupling: Taylor series
hamilt.sbc{1,1}.coeffs = 1;              % Linear coupling, slope 1

% Calculate (bound) eigen states (==> qm_bound)
hamilt.eigen.start = 00;                 % Lower index
hamilt.eigen.stop  = 20;                 % Upper index

% Save (bound) eigen states (==> qm_bound)
state.save_export = true;                 % Toggle saving to file(s)

% Initial density
time.rho.choice = 'thermal';             % thermal = Boltzmann density
time.rho.temperature = 0.1;              % choice of temperature

% Temperature and system-bath coupling
control.lvne.order = 'df';               % ordering vectorized density matrices: diagonals first
control.lvne.temperature = 0.1;          % temperature in atomic units: 315,923.5 K
control.relax.rate =  3.73779924e-4;     % Gamma_{2->0} should be = 1e-3 [D/hbar] 
control.relax.lower = 0;                 % Lower state for reference transition
control.relax.upper = 2;                 % Upper state for reference transition
control.relax.model = 'fermi';           % Omega dependent (Andrianov&Saalfrank)

% Define output observables and choose control targets
control.observe.types = 'prj';           % choose populations of observables
control.observe.choices = {0:2:10 1:2:9 11:20}; 
control.observe.labels = {'left well' 'right well' 'delocalized'};
control.observe.targets = 1:3;          % choose control targets 

% Parameters for balancing transformation
reduce.balance.A_stable = 'ssu';         % SSU or EVS method for stabilizing A  
reduce.balance.A_shift = 1e-6;           % magnitude of eigenvalue shift
reduce.balance.A_split = 1;              % dimension of instable part of A
reduce.balance.BN_scale = 250;           % Scaling factor for control fields
reduce.balance.acf_couple = true;        % additional control field coupling A and rho/x
reduce.balance.method = 'iter';          % ITER or BICG solver for GLYAP
reduce.balance.transform = 'srbt';       % SRBT or MRMR balancing transform;

% Parameters for H2 model reduction (Benner/Breiten @ MPI Magdeburg)
reduce.H2model.A_stable = 'ssu';         % SSU or EVS method for stabilizing A  
reduce.H2model.A_shift = 1e-6;           % magnitude of eigenvalue shift
reduce.H2model.A_split = 1;              % dimension of instable part of A
reduce.H2model.BN_scale = 150;           % Scaling factor for control fields
reduce.H2model.acf_couple = true;        % additional control field coupling A and rho/x
reduce.H2model.method = 'iter';          % ITER or BICG solver for BIRKA
reduce.H2model.conv_tol = 1e-6;          % Convergence tolerance for BIRKA
reduce.H2model.max_iter = 100;           % Max number of iterations for BIRKA

% Parameters for calculation of H2 error
reduce.H2error.A_shift = 1e-2;           % magnitude of eigenvalue shift
reduce.H2error.BN_scale = 10;            % Scaling factor for control fields
reduce.H2error.method = 'iter';          % ITER or BICG solver for GSYLV

% Plot of densities
if strcmpi(class(state),'wave')
    plots.density        = vis.contour;  % Curve plot of wave function
    plots.density.expect = false;
elseif strcmpi(class(state),'rho')
    plots.density = vis.bar;             % Bar graph of "ket" vector
end

% Plot of expectation values as time series
if strcmpi(class(state),'wave')
    plots.expect       = vis.expect;     % Normal representation
    plots.expect.e_max = 2;              % Range of energy plot
elseif strcmpi(class(state),'rho')
    plots.expect = vis.in_out;           % Input (fields u(t)) and output (observables y(t))
    plots.expect.y_min = -0.05;          % Range of population plot
    plots.expect.y_max = +0.80;          % Range of population plot
end