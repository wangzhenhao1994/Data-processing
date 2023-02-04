% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots space time

prt.disp ( '***************************************************************' )
prt.disp ( 'Vibrationally state selective excitation of a Morse ' )
prt.disp ( 'oscillator resembling one of the OH bonds in water  ' )
prt.disp ( 'Laser pulse optimized for 5-photon transition: 0->5 ' )
prt.disp ( '                                                    ' )
prt.disp ( 'see M.V.Korolkov, G.K.Paramonov, and B. Schmidt     ' )
prt.disp ( 'J. Chem. Phys. 105(5), 1862-1879  (1996), Fig. 2a   ' )
prt.disp ( '***************************************************************' )

% In the paper we used a reduced mass of 1728.539, but here 1728.257
m_O = atomic.mass.O16;
m_H = atomic.mass.H1;

% Spatial discretization
space.dof{1} = dof.fft;                  % Using fft grid
space.dof{1}.mass  = m_O*m_H/(m_O+m_H);  % Reduced mass 
space.dof{1}.n_pts = 064;                % Number of grid points
space.dof{1}.x_min =  0.7;               % Lower bound of grid 
space.dof{1}.x_max =  4.0;               % Upper bound of grid

% Hamiltonian operator 
hamilt.truncate.e_min = -0.1;            % Lower truncation of energy
hamilt.truncate.e_max =  0.5;            % Upper truncation of energy

hamilt.pot{1,1}     = pot.morse;         % Harmonic oscillator
hamilt.pot{1,1}.d_e = 0.1994;            % Dissociation energy
hamilt.pot{1,1}.r_e = 1.821;             % Equilibrium length
hamilt.pot{1,1}.alf = 1.189;             % Range parameter

hamilt.dip{1}{1,1} = dip.mecke;          % Mecke dipole function
hamilt.dip{1}{1,1}.r_0 = 0.6/atomic.l.A; % Length parameter: 0.6 A
hamilt.dip{1}{1,1}.q_0 = 7.85/atomic.d.D*atomic.l.A; % Charge parameter: 7.85 D/A

% Absorbing boundary conditions
hamilt.nip{1}      = nip.power;          % Negative imaginary potential
hamilt.nip{1}.exp  = 4;                  % Exponent
hamilt.nip{1}.min = +1.0;                % Beginning of inner grid region
hamilt.nip{1}.max = +3.5;                % End of inner grid region

% Temporal discretization
time.steps.m_delta  = 10/atomic.t.fs;    % Size of time steps: 10 fs 
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 100;               % Index of final time step
time.steps.s_number = 100;               % Number of sub steps per time step
% time.steps.s_number      = 1000;         % Number of sub steps per time step

% Electric field pulses
time.pulse{1}       = efi.sin_2;         % Shape of envelope
time.pulse{1}.polar = 0.0;               % Polarization angle [rad]
time.pulse{1}.delay = 500/atomic.t.fs;   % Time delay of pulse center
time.pulse{1}.fwhm  = 500/atomic.t.fs;   % Full width at half maximum
time.pulse{1}.ampli = 328.5/atomic.F.MV_cm;% From MV/cm to E_H/(e*a_0)
time.pulse{1}.frequ = 3424.19/atomic.w.cm_1;% From cm-1 to E_H

% Initial wave function
time.dof{1}     = init.morse;            % Ground state of Morse oscillator
time.dof{1}.m_r = space.dof{1}.mass;     % data copied from grid
time.dof{1}.d_e = hamilt.pot{1,1}.d_e;   % data copied from Morse potential
time.dof{1}.r_e = hamilt.pot{1,1}.r_e;   % data copied from Morse potential
time.dof{1}.alf = hamilt.pot{1,1}.alf;   % data copied from Morse potential
time.dof{1}.n_q = 0;                     % ground state

% Plot densities
plots.density        = vis.contour;      % Contour plot of Wigner transform
plots.density.expect = false;
plots.density.range  = true;             % Ranges of contour plot
plots.density.x_min  = 1;
plots.density.x_max  = 3.5;
plots.density.y_min  = -25;
plots.density.y_max  = 25;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 0.15;               % Range of energy plot
