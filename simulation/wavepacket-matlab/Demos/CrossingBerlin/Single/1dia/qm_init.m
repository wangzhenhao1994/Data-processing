% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global hamilt plots space time

prt.disp ( '***************************************************************')
prt.disp ( 'Single crossing example: Horenko, Salzmann, Schmidt, Schuette  ')
prt.disp ( '      see: The Journal of Chemical Physics, 117, 11075 (2002)  ')
prt.disp ( '                                                               ')
prt.disp ( '            ( A/R   C   )                                      ')
prt.disp ( 'V_dia (R) = (           )                                      ')
prt.disp ( '            ( C    BR^2 )                                      ')
prt.disp ( '                                                               ')
prt.disp ( 'with default values: A=1, B=1, C=0.1                           ')
prt.disp ( '                                                               ')
prt.disp ( 'Diabatic representation                                        ')
prt.disp ( '                                                               ')
prt.disp ( '***************************************************************')

% Number of (coupled) Schrödinger equations and representation
hamilt.coupling.n_eqs      = 2;          % Two coupled channels
hamilt.coupling.represent  = 'dia';      % Choose diabatic representation
hamilt.coupling.ini_rep    = 'dia';      % Choose diabatic initial state
hamilt.coupling.ini_coeffs = [1 0];      % Initially only upper adiabat

% Spatial discretization
space.dof{1}       = dof.fft;            % Using FFT based DVR/FBR scheme
space.dof{1}.mass  = 10000;              % Mass of particle
space.dof{1}.n_pts = 512;                % Number of grid points
space.dof{1}.x_min = 0.10;               % Lower bound of grid 
space.dof{1}.x_max = 2.85;               % Upper bound of grid

% Temporal discretization
time.steps.m_start = 000;                % Index of initial time step
time.steps.m_stop  = 100;                % Index of final time step
time.steps.m_delta = 1.0;                % Size of time steps 
time.steps.s_number = 050;               % Number of sub steps per time step

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width =  0.05;               % Width 
time.dof{1}.pos_0 =  0.4;                % Center in position representation
time.dof{1}.mom_0 =  100;                % Center in momentum representation

% Hamiltonian operator 
hamilt.truncate.e_min =  0.0;            % Lower truncation of energy
hamilt.truncate.e_max = 15.0;            % Upper truncation of energy

for m = 1:2
    for n = m:2
        hamilt.pot{m,n} = pot.single;    % Single Crossing
    end
end

% Absorbing boundary conditions
for m=1:2                                % Same for both channels
    hamilt.nip{m}      = nip.power;      % Negative imaginary potential
    hamilt.nip{m}.exp  = 4;              % Exponent
    hamilt.nip{m}.min = +0.15;           % Beginning of inner grid region
    hamilt.nip{m}.max = +2.50;           % End of inner grid region
end

% Plot time evolution of density
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
plots.expect.e_min    = 0;               % Manual setting of range for energies
plots.expect.e_max    = 5;               % Manual setting of range for energies
plots.expect.errorbar = false;           % Error bars for pot./kin. energy
plots.expect.legends  = true;            % Manual setting of range for population
plots.expect.export   = true;            % Export as JPG file

