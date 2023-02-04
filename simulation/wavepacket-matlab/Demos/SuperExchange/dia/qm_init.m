% Trying to reproduce figure 5 of Craig Martens' paper on
% Surface Hopping without Momentum Jumps: 
% A Quantum-Trajectory-Based Approach to Nonadiabatic Dynamics

% Copyright (C) 2020 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (k)
global hamilt plots space state time

prt.disp ( '***************************************************************')
prt.disp ( 'Trying to reproduce figure 5 of Craig Martens'' paper titled   ')
prt.disp ( 'Surface Hopping without Momentum Jumps:                        ')
prt.disp ( 'A Quantum-Trajectory-Based Approach to Nonadiabatic Dynamics   ')
prt.disp ( '                                                               ')
prt.disp ( '    Journal of Physical Chemistry A, 123(5), 1110-28 (2019)    ')
prt.disp ( '    DOI:10.1021/acs.jpca.8b10487                               ')
prt.disp ( '                                                               ')
prt.disp ( '***************************************************************')

% Number of (coupled) Schrödinger equations and representation
hamilt.coupling.n_eqs      = 3;          % Two coupled channels
hamilt.coupling.represent  = 'dia';      % Choose representation

% Wavepacket initially prepared in the first (=lowest) diabatic state
hamilt.coupling.ini_rep    = 'dia';      % Representation of initial state
hamilt.coupling.ini_coeffs = [1 0 0];    % Initial coefficients

% Spatial discretization
space.dof{1}       = dof.fft;            % Using FFT based DVR/FBR scheme
space.dof{1}.mass  = 2000;               % Mass of particle
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -25.0;              % Lower bound of grid 
space.dof{1}.x_max = +25.0;              % Upper bound of grid

% Temporal discretization
time.steps.m_start = 000;                % Index of initial time step
time.steps.m_stop  = 100;                % Index of final time step
time.steps.m_delta = 100;                % Size of time steps 
time.steps.s_number = 005;               % Number of sub steps per time step

% Surface hopping details (if applicable)
if isa (state,'sht.generic')
    state.rescale = 1;                   % Rescale momenta after hopping
    state.sca_nac = 0;                   % Rescale along NAC coupling vectors
end

% Superexchange model from page 3601 in DOI:10.1021/ct5003835
% Note that states 1 and 3 are not coupled directly
% while both 1 and 3 are coupled to state 2.
hamilt.pot{2,2}        = pot.taylor;     % constant potential 
hamilt.pot{2,2}.vshift = 0.010;
hamilt.pot{3,3}        = pot.taylor;     % constant potential 
hamilt.pot{3,3}.vshift = 0.005;
hamilt.pot{1,2}        = pot.gauss;      % Gaussian-shaped coupling
hamilt.pot{1,2}.pos_0  = 0.0;
hamilt.pot{1,2}.width  = sqrt(1/2);
hamilt.pot{1,2}.height = 0.001;
hamilt.pot{2,3}        = pot.gauss;      % Gaussian-shaped coupling
hamilt.pot{2,3}.pos_0  = 0.0;
hamilt.pot{2,3}.width  = sqrt(1/2);
hamilt.pot{2,3}.height = 0.010;

% Absorbing boundary conditions
for m=1:3                                % Same for all three channels
    hamilt.nip{m}      = nip.power;      % Negative imaginary potential
    hamilt.nip{m}.exp  = 4;              % Exponent
    hamilt.nip{m}.min = -20;             % Beginning of inner grid region
    hamilt.nip{m}.max = +20;             % End of inner grid region
end

% Initial wave function: Width is 10 times the inverse of the momentum k
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.pos_0 = -10;                 % Center in position representation
time.dof{1}.mom_0 =  k;                  % Center in momentum representation
time.dof{1}.width =  10/k;               % Width 

% Visualize densities evolving in time
plots.density           = vis.curve;     % Contour plot of the Wigner function
plots.density.represent = 'dvr';         % Position or momentum space
plots.density.complex   = 'abs2';        % Real or Imag or Abs2
plots.density.energy    = true;          % Show also energy function
plots.density.range     = false;         % Manual setting of plotting range

% plots.density = vis.surface;
% plots.density.srf_view  = [150 50];      % az el
% plots.density.srf_light = [25 75];       % az el
% plots.density.col_map = 'default';      
% plots.density.energy = true;
% plots.density.range    = true;           % Manual setting of plotting range
% plots.density.x_min = 0.2;
% plots.density.x_max = 2.7;
% plots.density.y_min = -50;
% plots.density.y_max = 280;
% plots.density.z_min = 00;
% plots.density.z_max = 10;

% plots.density = vis.contour;
% plots.density.marginals = true;
% plots.density.energy = true;
% plots.density.expect = false;
% plots.density.range    = true;           % Manual setting of plotting range
% plots.density.x_min = 0.10;
% plots.density.x_max = 2.85;
% plots.density.y_min = -50;
% plots.density.y_max = 280;

% Plot expectation values
plots.expect          = vis.expect; 
% plots.expect.e_min    = 0;               % Manual setting of range for energies
% plots.expect.e_max    = 5;               % Manual setting of range for energies
plots.expect.errorbar = false;           % Error bars for pot./kin. energy
plots.expect.legends  = true;            % Manual setting of range for population
plots.expect.export   = true;            % Export as JPG file

