% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global atomic hamilt plots space time

prt.disp ( '***************************************************************' )
prt.disp ( 'Vibrationally state selective electronic          ' )
prt.disp ( 'excitation of OH radical by a sin^2 pulse         ' )
prt.disp ( '***************************************************************' )
prt.disp ( '                                              ' )
prt.disp ( 'Light pulse optimized to induce a transition  ' )
prt.disp ( 'from the electronic *and* vibrational         ' )
prt.disp ( 'molecular ground state X2Pi(v=0)              ' )
prt.disp ( 'to an electronically *and* vibrationally      ' )
prt.disp ( 'molecular excited state A2Sigma+ (v=1)        ' )
prt.disp ( '                                              ' )
prt.disp ( 'see M.V.Korolkov, G.K.Paramonov               ' )
prt.disp ( 'Phys. Rev. A 57(6), 4998  (1998), Fig. 2a+b   ' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs      = 2;
hamilt.coupling.represent  = 'adi';
hamilt.coupling.labels     = {'X ^2\Pi', 'A ^2\Sigma^+'};
hamilt.coupling.ini_coeffs = [1 0];      % Initially: ground state only
hamilt.coupling.ini_rep    = 'dia';

% Atomic mases
m_O = atomic.mass.O16;
m_H = atomic.mass.H1;

% Spatial discretization
space.dof{1}       = dof.fft;            % using FFT grid
space.dof{1}.mass  = m_O*m_H/(m_O+m_H);  % Reduced mass (OH radical)
space.dof{1}.n_pts = 064;                % Number of grid points
space.dof{1}.x_min =  0.5;               % Lower bound of grid 
space.dof{1}.x_max =  5.0;               % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 100;               % Index of final time step
time.steps.m_delta  = 1/atomic.t.fs;     % Size of time steps: 1 fs 
time.steps.s_number = 100;               % Number of sub steps per time step

% Electric field as sequence of pulses
time.pulse{1}       = efi.sin_2;         % Shape of envelope
time.pulse{1}.delay = 50/atomic.t.fs;    % Time delay of pulse center: 50 fs
time.pulse{1}.fwhm  = 50/atomic.t.fs;    % Full width at half maximum: 50 fs
time.pulse{1}.ampli = 142.44/atomic.F.MV_cm; % Field amplitude: 142.44 MV/cm
time.pulse{1}.frequ = 35270.76/atomic.w.cm_1;% Carrier frequency: 35270.76 cm-1

% Hamiltonian operator 
hamilt.truncate.e_min = -0.2;            % Lower truncation of energy
hamilt.truncate.e_max =  0.5;            % Upper truncation of energy

hamilt.pot{1,1}     = pot.interp;        % Potentials: Interpolate table
hamilt.pot{2,2}     = pot.interp;        % Potentials: Interpolate table

hamilt.dip{1}       = cell(2);
hamilt.dip{1}{1,2}  = dip.interp;        % Dipoles: Interpolate table

% Absorbing boundary conditions
for m=1:2
    hamilt.nip{m}     = nip.power;       % Negative imaginary potential
    hamilt.nip{m}.exp = 4;               % Exponent
    hamilt.nip{m}.min = +1.0;            % Beginning of inner grid region
    hamilt.nip{m}.max = +4.5;            % End of inner grid region
end

% Initial wave function
time.dof{1}        = init.gauss;         % Gaussian-shaped wavepacket
time.dof{1}.width  =  0.1331;            % Width 
time.dof{1}.pos_0  =  1.8601;            % Center in position representation
time.dof{1}.mom_0  =  0.0;               % Center in momentum representation

% Plot densities
plots.density        = vis.contour;      % Contour plot
plots.density.expect = false;
plots.density.energy = true;
plots.density.range  = true;             % manually set the ranges
plots.density.x_min  = 1;
plots.density.x_max  = 3.5;
plots.density.y_min  = -20;
plots.density.y_max  = 20;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = +0.05;
plots.expect.e_min = -0.18;
