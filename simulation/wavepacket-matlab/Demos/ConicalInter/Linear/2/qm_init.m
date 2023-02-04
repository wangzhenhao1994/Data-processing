% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init ()
global hamilt plots space time

prt.disp ( '***************************************' )
prt.disp ( 'Linear J-T effect: Conical intersection' )
prt.disp ( 'adiabatic representation               ' )
prt.disp ( '***************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs      = 2;
hamilt.coupling.represent  = 'adi';
hamilt.coupling.ini_rep    = 'dia';      % Initial WF in diabatic representation
hamilt.coupling.ini_coeffs = [0 1];      % Initially only repulsive diabatic state populated

% Spatial discretization
space.dof{1}       = dof.fft;            % using FFT grid
space.dof{1}.mass  = 1;                  % mass
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max =  10.0;              % Upper bound of grid

% Spatial discretization
space.dof{2}       = dof.fft;            % using FFT grid
space.dof{2}.mass  = 1;                  % mass
space.dof{2}.n_pts = 64;                 % Number of grid points
space.dof{2}.x_min = -10.0;              % Lower bound of grid 
space.dof{2}.x_max =  10.0;              % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 019;               % Index of final time step
time.steps.m_delta  = 0.10;              % Size of time steps 
time.steps.s_number = 0050;              % Number of sub steps per time step

% Hamiltonian operator 
hamilt.truncate.e_min    = -15.0;        % Lower truncation of energy
hamilt.truncate.e_max    =  50.0;        % Upper truncation of energy

for m = 1:hamilt.coupling.n_eqs
    for n = m:hamilt.coupling.n_eqs
        hamilt.pot{m,n} = pot.con_int;   % Conical intersection
        hamilt.pot{m,n}.omega=  0.0;     % Harmonic frequency
        hamilt.pot{m,n}.kappa=  1.0;     % Linear JT coupling
        hamilt.pot{m,n}.gamma=  0.0;     % Quadratic JT coupling
    end
end

for m = 1:hamilt.coupling.n_eqs
    hamilt.nip{m} = nip.power;           % Negative imaginary potential
    hamilt.nip{m}.exp  = [4 4];          % Exponent
    hamilt.nip{m}.min  = [-08.0 -08.0];  % Lower bound
    hamilt.nip{m}.max  = [+08.0 +08.0];  % Upper bound
end

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width =  0.7;                % Width 
time.dof{1}.pos_0 = -4.0;                % Center in position representation
time.dof{1}.mom_0 =  4.0;                % Center in momentum representation

time.dof{2}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{2}.width =  1.0;                % Width 
time.dof{2}.pos_0 =  0.0;                % Center in position representation
time.dof{2}.mom_0 =  0.0;                % Center in momentum representation

% Plot densities
plots.density           = vis.surface;
plots.density.represent = 'dvr';         % Position (dvr) or momentum (fbr) densities
plots.density.energy    = true;          % Show also energy functions
plots.density.srf_view  = [15 42];       % View point for surface plot: [az el]
plots.density.srf_look  = [false false]; % Look of surface plot [shading lighting]
plots.density.srf_light = [45 45];       % Angle of light for surface plot: [az el]

% plots.density           = vis.contour;
% plots.density.range = true;              % manual setting of plotting range
% plots.density.x_min = -7;
% plots.density.x_max = 8.5;
% plots.density.y_min = -4;
% plots.density.y_max = 4;

% Plot expectation values
plots.expect            = vis.expect;
plots.expect.e_max      = 20;            %
