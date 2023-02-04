% Copyright (C) 2018-20.. Burkhard Schmidt's group
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_test (n_dim, n_per)
global hamilt space time % plots

prt.disp ( '***************************************************************' )
prt.disp ( 'Benchmark: Coherent state of a harmonic oscillator, N dimensions' )
prt.disp ( '***************************************************************' )

% Spatial discretization
for i_dim=1:n_dim
    space.dof{i_dim} = dof.fft;          % Using FFT grid method
    space.dof{i_dim}.n_pts = 032;        % Number of grid points
    space.dof{i_dim}.x_min = -10.0;      % Lower bound of grid
    space.dof{i_dim}.x_max =  10.0;      % Upper bound of grid
    space.dof{i_dim}.mass = 1;           % Mass
end

% Temporal discretization
time.steps.m_start = 000;                % 
time.steps.m_stop  = 10*n_per;           % FIVE PERIODS OF OSCILLATION
time.steps.m_delta = pi/05;              % Size of time steps 
time.steps.s_number = 100;               % Number of sub steps per time step

% Initial state
for i_dim=1:n_dim
    time.dof{i_dim}       = init.gauss;  % Gaussian-shaped wavepacket
    time.dof{i_dim}.width = sqrt(1/2);   % Width
    time.dof{i_dim}.pos_0 = -5.0;        % Center in position representation
    time.dof{i_dim}.mom_0 =  0.0;        % Center in momentum representation
end

% Truncate Hamiltonian operator 
hamilt.truncate.e_min  =  0.0;           % Lower truncation of energy
hamilt.truncate.e_max  = 50 * n_dim;     % Upper truncation of energy

% Potential energy
c1 = zeros(1,n_dim);
c2 = ones(1,n_dim);
hamilt.pot{1,1}        = pot.taylor;     % Taylor expansion
hamilt.pot{1,1}.coeffs = [c1;c2];        % Force constant

% Plot time evolution of density
% plots.density            = vis.contour;
% plots.density.marginals  = true;
% plots.density.energy     = true;
% plots.density.expect     = true;

% Plot expectation values
% plots.expect          = vis.expect;
% plots.expect.errorbar = false;           % Error bars for pot./kin. energy
% plots.expect.legends  = true;            % Manual setting of range for population
% plots.expect.export   = true;            % Export as JPG file

% Plot spectrum
% plots.spectrum = vis.spectrum;          % Power spectrum
% plots.spectrum.export = true;            % Export as JPG file




