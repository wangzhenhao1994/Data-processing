% Copyright (C) 2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~)
global atomic hamilt plots space time

prt.disp ('********************************')
prt.disp ('Selective bond breakage of HOD  ')
prt.disp ('with a shaped laser pulse.      ')
prt.disp ('see PhysRev A 78:065402 (2008)  ')
prt.disp ('********************************')

% Number of coupled equations
hamilt.coupling.n_eqs      = 2;
hamilt.coupling.represent  = 'dia';
hamilt.coupling.labels     = {'X^1 A_1', 'A^1 B_1'};
hamilt.coupling.ini_coeffs = [1 0];      % Population of initial state
hamilt.coupling.ini_rep    = 'dia';      % Representation of initial wave function

% Grid definition; masses are not needed here since we use
% a custom operator; it would be smarter to use an unbalanced
% grid (O-D is heavier and thus needs more grid points than O-H),
% but we shall stick to the original setup as close as possible.
space.dof{1}       = dof.fft;            % using FFT grid
space.dof{1}.label = 'O-D';              % label for axis descriptions
space.dof{1}.n_pts = 512;                % number of points
space.dof{1}.x_min = 1.15;               % lower boundary of the grid
space.dof{1}.x_max = 26.75;              % upper boundary of the grid

space.dof{2}       = dof.fft;            % using FFT grid
space.dof{2}.label = 'O-H';              % label for axis descriptions
space.dof{2}.n_pts = 512;                % number of points
space.dof{2}.x_min = 1.15;               % lower boundary of the grid
space.dof{2}.x_max = 26.75;              % upper boundary of the grid

hamilt.amo{1} = amo.reaction;            % projecting on D + OH channel
hamilt.amo{1}.reac  = 1;
hamilt.amo{1}.prod  = 2;
hamilt.amo{1}.side  = 'r';
hamilt.amo{1}.label = 'D + OH';

hamilt.amo{2} = amo.reaction;            % projecting on H + OD channel
hamilt.amo{2}.reac  = 1;
hamilt.amo{2}.prod  = 2;
hamilt.amo{2}.side  = 'p';
hamilt.amo{2}.label = 'H + OD';

% Temporal discretisation
time.steps.m_start  = 0;                 % index of first time step
time.steps.m_stop   = 100;               % index of last time step
time.steps.m_delta  = 20;                % approx 0.5 fs per time step
time.steps.s_number 

% Electric field loaded from file
time.pulse{1}         = efi.inter;       % Interpolate from file; We interpolate the
time.pulse{1}.ampli   = 1;               % whole pulse from file, including the
time.pulse{1}.delay   = 0;               % oscillations, so all other variables only get
time.pulse{1}.frequ   = 0;               % dummy values.
time.pulse{1}.file    = 'field_x.dat';   % file to load the field from
time.pulse{1}.method  = 'spline';        % Spline interpolation

% Hamiltonian operator
hamilt.truncate.e_min = -0.5;            % lower truncation of energy
hamilt.truncate.e_max = 0.5;             % upper truncation of energy

hamilt.kin{1}       = kin.triatomic;     % Kinetic energy for fixed bending angle
hamilt.kin{1}.theta = 104.52 * pi/180;   % bending angle
hamilt.kin{1}.dof   = [1 2];             % indices of coordinates for AB, BC distance
hamilt.kin{1}.mass  = [2.13 16.00 1.00] / atomic.m.u;
                                         % note that the D mass is off by a few percent
                                         % due to approximations in Ashwanis setup
                                         % (he assumed mu_OD = 2 * mu_OH).

% Potential energy functions
hamilt.pot{1,1} = pot.H2O;               % potential energy (ground state)
hamilt.pot{2,2} = pot.H2O;               % potential energy (excited state)

% Transition dipole moment
hamilt.dip{1}      = cell(2);
hamilt.dip{1}{1,2} = dip.H2O;            % transition dipole moments

% Absorbing boundary conditions
hamilt.nip{1} = nip.power;               % Negative imaginary potential
hamilt.nip{1}.exp = [2 2];               % Exponent
hamilt.nip{1}.min = [1.3 1.3];           % Begin of inner grid region
hamilt.nip{1}.max = [25 25];             % End of inner grid region

hamilt.nip{2} = nip.power;               % Negative imaginary potential
hamilt.nip{2}.exp = [2 2];               % Exponent
hamilt.nip{2}.min = [1.3 1.3];           % Begin of inner grid region
hamilt.nip{2}.max = [25 25];             % End of inner grid region

% Initial wave function; HOD ground state modelled by two Gaussians
time.dof{1}        = init.gauss;         % Gaussian initial wave function
time.dof{1}.width  = 1/sqrt(4*23.18);    % width of Gaussian
time.dof{1}.pos_0  = 1.806;              % center in position representation
time.dof{1}.mom_0  = 0;                  % center in momentum representation

time.dof{2} = init.gauss;                % Gaussian initial wave function
time.dof{2}.width  = 1/sqrt(4*18.84);    % Width of Gaussian
time.dof{2}.pos_0  = 1.833;              % Center in position representation
time.dof{2}.mom_0  = 0;                  % center in momentum representation

% Plot time evolution of density
plots.density        = vis.contour;      % Draw contour lines
plots.density.expect = false;
plots.density.range  = true;             % manually set ranges
plots.density.x_min  = 1.2;
plots.density.x_max  = 6;
plots.density.y_min  = 1.2;
plots.density.y_max  = 6;

plots.density.cnt_nlev = [50 15];        % number of contour lines
plots.density.scale_dvr = 6;             % factor that determines ranges of density
                                         % for contour lines to draw
plots.density.pot_min = -0.35;           % finetuning for drawing of energy curves
plots.density.pot_max = 0;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 0.3;                % manually set ranges of energy plot