% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic control hamilt plots space state time

prt.disp ( '***************************************************************' )
prt.disp ( 'Vibrationally state selective excitation of a Morse ' )
prt.disp ( 'oscillator resembling one of the OH bonds in water  ' )
prt.disp ( 'Decay of 5th vibrationally excited state            ' )
prt.disp ( '                                                    ' )
prt.disp ( 'see M.V.Korolkov, G.K.Paramonov, and B. Schmidt     ' )
prt.disp ( 'J. Chem. Phys. 105(5), 1862-1879  (1996)            ' )
prt.disp ( '***************************************************************' )

state.save_export = true; % Save wavefunctions etc

% In the paper we used a reduced mass of 1728.539, but here 1728.257
m_O = atomic.mass.O16;
m_H = atomic.mass.H1;

% Spatial discretization
space.dof{1} = dof.fft;                  % Using fft grid
space.dof{1}.mass  = m_O*m_H/(m_O+m_H);  % Reduced mass 
space.dof{1}.n_pts = 512;                % Number of grid points
space.dof{1}.x_min = 0.7;                % Lower bound of grid 
space.dof{1}.x_max = 12.0;               % Upper bound of grid

% Hamiltonian operator 
hamilt.truncate.e_min = -0.1;            % Lower truncation of energy
hamilt.truncate.e_max =  0.5;            % Upper truncation of energy

hamilt.pot{1,1}      = pot.morse;        % Harmonic oscillator
hamilt.pot{1,1}.d_e  = 0.1994;           % Dissociation energy
hamilt.pot{1,1}.r_e  = 1.821;            % Equilibrium length
hamilt.pot{1,1}.alf  = 1.189;            % Range parameter

hamilt.dip{1}{1,1} = dip.mecke;          % Mecke dipole function
hamilt.dip{1}{1,1}.r_0 = 0.6/atomic.l.A; % Length parameter: 0.6 A
hamilt.dip{1}{1,1}.q_0 = 7.85/atomic.d.D*atomic.l.A; % Charge parameter: 7.85 D/A

hamilt.sbc{1,1} = sbc.taylor;            % System-bath coupling: Taylor series
hamilt.sbc{1,1}.coeffs = 1;              % Linear coupling, slope 1

% Temporal discretization
time.steps.m_delta  = 10/atomic.t.fs;    % Size of time steps: 10 fs 
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 100;               % Index of final time step
time.steps.s_number = 100;               % Number of sub steps per time step

% Calculate and save (bound) eigen states (==> qm_bound)
hamilt.eigen.stop = 00;                  % Lower index
hamilt.eigen.stop = 21;                  % Upper index

% Define populations as observables
control.observe.types = 'prj';           % types of observables
control.observe.choices = {[0] [1] [2] [3] [4] [5]};
control.observe.labels  = {'|0>' '|1>' '|2>' '|3>' '|4>' '|5>'};
control.observe.targets = 1:6;           % evaluate targets at terminal time

% Initial state (==> qm_abncd)
time.rho.choice = 'pure';                % starting from a pure state
time.rho.pure = 5;                       % 5-th vibrationally excited state

% Plot densities: contours or bar
if strcmpi(class(state),'wave')
    plots.density = vis.curve;           % Curve plot of wave function
elseif strcmpi(class(state),'rho')
    plots.density = vis.bar;             % Bar graph of "ket" vector
end

% Plot expectation values
if strcmpi(class(state),'wave')
    plots.expect = vis.expect;           % "normal" representation
    plots.expect.e_min = -0.05;          % Range of energy plot
    plots.expect.e_max = +0.25;          % Range of energy plot
elseif strcmpi(class(state),'rho')
    plots.expect = vis.in_out;           % input (fields u(t)) and output (observables y(t))
    plots.expect.y_min = -0.10;          % Range of population plot
    plots.expect.y_max = +1.10;          % Range of population plot
end

% Open quantum system (LvNE): temperature, relaxation, (dephasing)
control.lvne.temperature = 0.00;         % Temperature in atomic units: 315,923.5 K
control.lvne.order = 'df';               % Split into densities and coherences
control.relax.model = 'fermi';           % Choice of relaxation model
control.relax.rate  = 2*atomic.t.ps;     % Relaxation rate (inverse of 1/e time)
control.relax.lower = 0;                 % Lower state for reference transition
control.relax.upper = 1;                 % Upper state for reference transition
% control.depha.model = 'gauss';           % Choice of relaxation model
% control.depha.rate  = 2*atomic.t.ps;     % Relaxation rate (inverse of 1/e time)
% control.depha.lower = 0;                 % Lower state for reference transition
% control.depha.upper = 1;                 % Upper state for reference transition