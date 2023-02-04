% Copyright (C) 2004-201Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global hamilt plots space state time

prt.disp ( '***************************************************************')
prt.disp ( 'Dual crossing example :  Horenko, Salzmann, Schmidt, Schuette  ')
prt.disp ( '      see: The Journal of Chemical Physics, 117, 11075 (2002)  ')
prt.disp ( '                                                               ')
prt.disp ( '             ( A     C   )                                     ')
prt.disp ( ' V_dia (R) = (           )                                     ')
prt.disp ( '             ( C    BR^2 )                                     ')
prt.disp ( '                                                               ')
prt.disp ( 'with default values: A=1, B=1, C=0.1                           ')
prt.disp ( '                                                               ')
prt.disp ( 'Adiabatic representation                                       ')
prt.disp ( '                                                               ')
prt.disp ( '***************************************************************')

% Number of (coupled) Schrödinger equations and representation
hamilt.coupling.n_eqs      = 2;
hamilt.coupling.represent  = 'adi';
hamilt.coupling.ini_rep    = 'adi';
hamilt.coupling.ini_coeffs = [1 0];      % Initially only lower adiabatic state populated

% Spatial discretization
space.dof{1}      = dof.fft;             % using FFT grid
space.dof{1}.mass = 2000;                % mass
space.dof{1}.n_pts = 0512;               % Number of grid points
space.dof{1}.x_min = -4.0;               % Lower bound of grid 
space.dof{1}.x_max =  5.0;               % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 100;               % Index of final time step
time.steps.m_delta  = 5.0;               % Size of time steps 
time.steps.s_number = 050;               % Number of sub steps per time step

% Surface hopping details (if applicable)
if isa (state,'sht.generic')
    state.rescale = 1;                   % Rescale momenta after hopping
    state.sca_nac = 0;                   % Rescale along NAC coupling vectors
end

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width =  0.2991/2;           % Width 
time.dof{1}.pos_0 = -2.5;                % Center in position representation
time.dof{1}.mom_0 = 20.0;                % Center in momentum representation

% Hamiltonian operator 
hamilt.truncate.e_min    = -0.5;         % Lower truncation of energy
hamilt.truncate.e_max    = +5.0;         % Upper truncation of energy

% Potential energy
% hamilt.pot{1,1}        = pot.taylor;     % Taylor series
% hamilt.pot{1,1}.vshift = 1;              % Constant potential
% hamilt.pot{2,2}        = pot.taylor;     % Taylor series
% hamilt.pot{2,2}.coeffs = [0;2];          % Parabolic potential
% hamilt.pot{1,2}        = pot.taylor;     % Taylor series
% hamilt.pot{1,2}.vshift = 0.1;            % Coupling (constant)

for m = 1:2
    for n = m:2
        hamilt.pot{m,n} = pot.dual;      % Dual crossing
    end
end

% Absorbing boundary conditions
for m = 1:2
    hamilt.nip{m}      = nip.power;      % Negative imaginary potential
    hamilt.nip{m}.exp  = 4;              % Exponent
    hamilt.nip{m}.min = -4.0;            % Beginning of inner grid region
    hamilt.nip{m}.max = +4.0;            % End of inner grid region
end

% Plot densities
plots.density           = vis.curve;     % Contour plot of the Wigner function
plots.density.represent = 'dvr';         % Position or momentum space
plots.density.complex   = 'abs2';        % Real or Imag or Abs2
plots.density.energy    = true;          % Show also energy function
plots.density.range     = false;         % Manual setting of plotting range

% plots.density        = vis.contour;
% plots.density.expect = false;
% plots.density.range  = true;             % manual setting of plotting range
% plots.density.x_min  = -3.0;
% plots.density.x_max  = +4.5;
% plots.density.y_min  = -30;
% plots.density.y_max  = +80;
% % plots.density.srf_view = [-7 35];        % View point for surface plot: [az el]

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 1.5;
