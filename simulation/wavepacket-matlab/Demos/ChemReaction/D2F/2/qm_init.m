% Copyright (C) 2008-2009 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global atomic hamilt plots space time

prt.disp ( '********************************************************************' )
prt.disp ( 'Collinear model for the reaction F+D2 -> DF + D        ' )
prt.disp ( 'see also Whitlock, Muckerman, JCP 61:4618 (1974) and   ' )
prt.disp ( 'Muckerman, JCP 54:1155 (1971)                          ' )
prt.disp ( '                                                       ' )
prt.disp ( 'LEPS surface taken from Steinfeld, Francisco, Hase     ' )
prt.disp ( '"Chemical Kinetics and Dynamics" p. 245, problem 7.14  ' )
prt.disp ( '                                                       ' )
prt.disp ( 'High collision energy                                  ' )
prt.disp ( '********************************************************************' )

% atomic masses
m_F = atomic.mass.F19;                   % mass of fluorine (a.u.)
m_D = atomic.mass.H2;                    % mass of deuterium(a.u.)

% Spatial discretization
space.dof{1}       = dof.fft;            % D-D distance, equally spaced grid
space.dof{1}.label = 'D-D';
space.dof{1}.mass  = 1;                  % the internal kinetic energy is disabled
space.dof{1}.x_min = 0.4;                % minimum D-D distance
space.dof{1}.x_max = 7;                  % maximum D-D distance
space.dof{1}.n_pts = 128;                % number of points

space.dof{2}       = dof.fft;            % D-F distance, equally spaced grid
space.dof{2}.label = 'D-F';
space.dof{2}.mass  = 1;                  % internal KE is disabled anyway
space.dof{2}.x_min = 0.6;                % minimum D-F distance
space.dof{2}.x_max = 8;                  % maximum D-F distance
space.dof{2}.n_pts = 128;                % number of points

hamilt.amo{1} = amo.reaction;            % population of reactant states
hamilt.amo{1}.reac = 2;                  % educt distance is DOF #2
hamilt.amo{1}.prod = 1;                  % product distance is DOf #1
hamilt.amo{1}.side = 'r';                % projection on the reactants
hamilt.amo{1}.label= 'D_2 + F';          % labeling the reactants

hamilt.amo{2} = amo.reaction;            % population of product states
hamilt.amo{2}.reac = 2;                  % educt distance is DOF #2
hamilt.amo{2}.prod = 1;                  % product distance is DOf #1
hamilt.amo{2}.side = 'p';                % projection on the products
hamilt.amo{2}.label= 'DF + D';           % labeling the products

% Temporal discretization
time.steps.m_start  = 0;                 % first time step
time.steps.m_stop   = 100;               % last time step
time.steps.m_delta  = 1/atomic.t.fs;     % length of one time step
time.steps.s_number = 50;                % number of sub steps

% Hamiltonian operator 
hamilt.truncate.e_min= -0.3;             % lower truncation of energy
hamilt.truncate.e_max=  0.3;             % upper truncation of energy

hamilt.kin{1}      = kin.triatomic;      % triatomic kinetic energy operator
hamilt.kin{1}.dof  = [1 2];              % act on the first two coordinates
hamilt.kin{1}.mass = [m_D m_D m_F];      % masses of the particles
hamilt.kin{1}.theta = pi;                % bending angle (pi for linear case)

hamilt.pot{1,1}     = pot.leps;               % LEPS surface
hamilt.pot{1,1}.ang = pi;                     % bending angle
hamilt.pot{1,1}.d_e = [0.1745 0.2251 0.2251]; % diss. energy
hamilt.pot{1,1}.r_e = [1.402  1.7329 1.7329]; % equilib. dist.
hamilt.pot{1,1}.alf = [1.0277 1.1742 1.1742]; % range parameter
hamilt.pot{1,1}.s_p = [0.106  0.167  0.167 ]; % Sato parameter

% Absorbing boundary conditions
hamilt.nip{1}  = nip.power;              % Absorbing boundary conditions
hamilt.nip{1}.exp  = [2 2];
hamilt.nip{1}.min  = [0.6 0.8];
hamilt.nip{1}.max  = [6 7];

% Initial wave function
% We assume D-D in the HO ground state, and F coming in at a
% low speed (E_kin ~ 25 meV).
%
% Note that the coordinate transformation also affects the
% initial wave function!
%
% In atomic positions (x_a = Fluorine, x_b = first deuterium, x_c = second deuterium)
% the initial wave function is simply
% 
% exp(i k_0 x_a) * exp(-i k_0 0.5*(x_b+x_c) * exp(-a(1/2(x_b-x_c) - R0)^2)
%                * f(x_a) * g(0.5(x_b + x_c))
%
% Read: D2 and F crash into each other with momentum k_0 each (CMS frame), 
% D2 is in the vibrational ground state, and F and D2 are located somewhere
% (described by the functions f and g). If we transform this to bond length
% coordinates, we find that the first two terms give
%
% exp(-i k_0 R_ab) * exp(-i k_0/2 R_bc),
%
% i.e., the initial wave function for the D-D distance needs to get half
% the momentum, too. The transformations of f and g are nontrivial, too.
% Here, we just assume a Gaussian for the D-F distance.

omega = sqrt(2*0.174*1.0277^2/(m_D/2));

time.dof{1}        = init.gauss;
time.dof{1}.pos_0  = 1.4019;
time.dof{1}.mom_0  = -3.5;
time.dof{1}.width  = sqrt(1/(2*(m_D/2)*omega));

time.dof{2}        = init.gauss;
time.dof{2}.pos_0  = 5;
time.dof{2}.mom_0  = -7;
time.dof{2}.width  = 0.3;

% Plot time evolution of density
plots.density          = vis.surface; 
plots.density.srf_view = [110 70];

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 0.2;
