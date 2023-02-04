%------------------------------------------------------------------------------
%
% Create an object "state" before starting a new WavePacket calculation
% where for "state" there is a substantial choice of different classes.
% Also clears the workspace to avoid side-effects from multiple runs.
% Also closes/recreates log files etc. Call this function whenever 
% you need a fresh workspace.
%
% This function also provides numeric conversion factors between
% atomic units (used inside WavePacket) and SI (and related) units
%
% https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_setup
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-20.. Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2008 Burkhard Schmidt
%               2011 Ulf Lorenz, Boris Schaefer-Bung
%               2012 Jeremy Rodriguez, Ulf Lorenz
%               2016 Burkhard Schmidt
%
% see the README file for license details.

function qm_setup (string1, param1, param2, param3)

% Clear variables/functions and close files (from previous runs);
clear global atomic correlate control expect hamilt info plots space state time
% clear functions;
% clear classes;
fclose('all');

% Using global variables
global info state time

% Running under Matlab or Octave
if exist ('OCTAVE_VERSION', 'builtin')
    info.system = 'Octave';
    warning("off", "Octave:shadowed-function");
    pkg load statistics
else
    info.system = 'Matlab';
end


% If we recycle plots under Linux Octave, they may sometimes be screwed up somehow.
% Simply close all of them.
if strcmpi(info.system, 'Octave')
    close all
end


% Choice of WavePacket main classes 
if nargin==0
    string1 = 'wave';                    % Default: wavefunctions
end

switch lower(string1)
    case 'wave'     % quantum wavefunctions, represented on grids    
        state = wave;
        
    case 'traj'     % phase-space densities, sampled by classical trajectories
        if nargin==1
            param1=[];                   % number of trajectories
            param2=[];                   % seeding random number generator
        elseif nargin==2
            param2=[];
        end
        state = traj(param1,param2); 
        
    case 'ket'      % quantum state vectors in eigen|energy representation
        state = ket;
        
    case 'rho'      % quantum density matrix in eigen|energy representation
        state = rho;
    
    case 'fssh'     % fewest switches surface hopping
        if nargin==1
            param1=[];                   % number of trajectories
            param2=[];                   % seeding random number generator
        elseif nargin==2
            param2=[];
        end
        state = sht.fssh(param1,param2);
        
    case 'qtsh'     % fewest switches surface hopping
        if nargin==1
            param1=[];                   % number of trajectories
            param2=[];                   % seeding random number generator
        elseif nargin==2
            param2=[];
        end
        state = sht.qtsh(param1,param2);
     
    case 'mssh'     % multiple switches surface hopping
        if nargin==1
            param1=[];                   % number of trajectories
            param2=[];                   % seeding random number generator
        elseif nargin==2
            param2=[];
        end
        state = sht.mssh(param1,param2);
    
    case {'sssh','lz'}                   % single switch surface hopping
        if nargin==1
            param1=[];                   % LZ variant
            param2=[];                   % number of trajectories
            param3=[];                   % seeding random number generator
        elseif nargin==2
            param2=[];                   % number of trajectories
            param3=[];                   % seeding random number generator
        elseif nargin==3
            param3=[];                   % seeding random number generator
        end
        state = sht.sssh(param1,param2,param3);
        
    case {'sssh_gd','lz_gd'}             % single switch surface hopping, with gradient descent
        if nargin==1
            param1=[];                   % number of trajectories
            param2=[];                   % seeding random number generator
        elseif nargin==2
            param2=[];                   % seeding random number generator
        end
        if nargin==1
            param1=[];                   % LZ variant
            param2=[];                   % number of trajectories
            param3=[];                   % seeding random number generator
        elseif nargin==2
            param2=[];                   % number of trajectories
            param3=[];                   % seeding random number generator
        elseif nargin==3
            param3=[];                   % seeding random number generator
        end
        state = sht.sssh_gd(param1,param2,param3);
        
    case {'sssh_gd_mc','lz_gd_mc'}        % what does "multiple channels" exactly mean?
        if nargin==1
            param1=[];                   % number of trajectories
            param2=[];                   % seeding random number generator
        elseif nargin==2
            param2=[];                   % seeding random number generator
        end
        if nargin==1
            param1=[];                   % LZ variant
            param2=[];                   % number of trajectories
            param3=[];                   % seeding random number generator
        elseif nargin==2
            param2=[];                   % number of trajectories
            param3=[];                   % seeding random number generator
        elseif nargin==3
            param3=[];                   % seeding random number generator
        end
         state = sht.sssh_gd_mc(param1,param2,param3);
        
    otherwise
        prt.error ('Wrong choice of WavePacket main class')
end

% Display (recently constructed) object "state"
disp(state)

% Calling a few constructors: temporal discretization
time.steps  = tmp.steps;
time.efield = tmp.efield;

% Calling a few constructors: setting up the Hamiltonian
global hamilt
hamilt.coupling = ham.coupling;
hamilt.eigen    = ham.eigen;
hamilt.truncate = ham.truncate;

%% Converting between SI units and atomic units
global atomic

% Constants in SI units (taken from Wikipedia, 2016)
m_e = 9.10938215E-31;                    % electron rest mass [kg]
q_e = 1.60217662E-19;                    % elementary charge [C]
h_b = 1.05457186E-34;                    % Planck's constant [Js] 
k_e = 8.9875517873682E+9;                % electric constant [(kg*m^3)/(s^2*C^2)]
n_a = 6.0221408E23;                      % Avogadro's number [1/mol]
c_v = 299792458;                         % vacuum speed of light [m/s]
k_b = 1.380648E-23;                      % Boltzmann's constant [J/K]

% Atomic unit of length
atomic.m.kg = m_e;                       % rest mass of electron [kg]
atomic.m.g  = m_e * 1e3;                 % same, but in gram 
atomic.m.u  = atomic.m.g * n_a;          % same, but in g/mol 

% Atomic unit of length
atomic.l.m = h_b^2 / (k_e*m_e*q_e^2);    % Bohr radius [m]
atomic.l.nm = atomic.l.m * 1e09;         % nanometer
atomic.l.A =  atomic.l.m * 1e10;         % Angstrom
atomic.l.pm = atomic.l.m * 1e12;         % picometer

% Atomic unit of energy (and equivalents)
atomic.E.J = m_e*q_e^4*k_e^2 / h_b^2;    % Hartree energy [J]
atomic.E.J_M  = atomic.E.J * n_a;        % Joule per mole
atomic.E.kJ_M = atomic.E.J_M / 1e3;      % kiloJoule per mole    
atomic.E.c_M  = atomic.E.J_M / 4.184;    % Calorie per mole
atomic.E.kc_M = atomic.E.c_M / 1e3;      % kiloCalorie per mole
atomic.E.eV   = atomic.E.J / q_e;        % ElectronVolt
atomic.E.meV  = atomic.E.eV * 1e3;       % MilliElectronVolt
atomic.E.ueV  = atomic.E.eV * 1e6;       % MicroElectronVolt

% Atomic unit of temperature
atomic.T.K  = atomic.E.J / k_b;          % Kelvin
atomic.T.mK = atomic.T.K * 1e3;          % milliKelvin
atomic.T.uK = atomic.T.K * 1e6;          % microKelvin

% Atomic unit of wavenumbers
atomic.w.m_1  = atomic.E.J / (2*pi*h_b*c_v);  % Inverse Meter  
atomic.w.cm_1 = atomic.w.m_1 / 100;      % Inverse CentiMeter    

% Atomic unit of time
atomic.t.s  = h_b / atomic.E.J;          % Second
atomic.t.ms = atomic.t.s * 1e03;         % MilliSecond
atomic.t.us = atomic.t.s * 1e06;         % MicroSecond
atomic.t.ns = atomic.t.s * 1e09;         % NanoSecond
atomic.t.ps = atomic.t.s * 1e12;         % PicoSecond
atomic.t.fs = atomic.t.s * 1e15;         % FemtoSecond
atomic.t.as = atomic.t.s * 1e18;         % AttoSecond

% Atomic unit of (angular!) frequency "omega"
atomic.o.Hz =  atomic.E.J / h_b;         % Hertz
atomic.o.kHz = atomic.o.Hz / 1e03;       % KiloHertz
atomic.o.MHz = atomic.o.Hz / 1e06;       % MegaHertz
atomic.o.GHz = atomic.o.Hz / 1e09;       % GigaHertz
atomic.o.THz = atomic.o.Hz / 1e12;       % TeraHertz
atomic.o.PHz = atomic.o.Hz / 1e15;       % PetaHertz

% Atomic unit of dipole moment
atomic.d.Cm = atomic.l.m * q_e;          % Coulomb * Meter 
atomic.d.D  = atomic.d.Cm * c_v * 1e21;  % Debye

% Atomic unit of electric field
atomic.F.V_m   = atomic.E.J / (q_e * atomic.l.m); 
atomic.F.MV_m  = atomic.F.V_m / 1e06;    % MegaVolt per Meter
atomic.F.MV_cm = atomic.F.V_m / 1e08;    % MegaVolt per CentiMeter
atomic.F.GV_m  = atomic.F.V_m / 1e09;    % GigaVolt per Meter
atomic.F.GV_cm = atomic.F.V_m / 1e11;    % GigaVolt per CentiMeter
atomic.F.V_A   = atomic.F.V_m / 1e10;    % Volt per Angstrom

% Atomic unit of light intensity
atomic.I.W_m2   = atomic.F.V_m^2*c_v/(8*pi*k_e);  % Watt per Meter^2
atomic.I.W_cm2  = atomic.I.W_m2  / 1e04; %     Watt per CentiMeter^2
atomic.I.kW_cm2 = atomic.I.W_cm2 / 1e03; % KiloWatt per CentiMeter^2
atomic.I.MW_cm2 = atomic.I.W_cm2 / 1e06; % MegaWatt per CentiMeter^2
atomic.I.GW_cm2 = atomic.I.W_cm2 / 1e09; % GigaWatt per CentiMeter^2
atomic.I.TW_cm2 = atomic.I.W_cm2 / 1e12; % TeraWatt per CentiMeter^2
atomic.I.PW_cm2 = atomic.I.W_cm2 / 1e15; % PetaWatt per CentiMeter^2

%% Table of atomic masses (from Wikipedia, webelements.com)
atomic.mass.neutron  = 1.008664916/atomic.m.u; 
atomic.mass.electron = 5.485799091E-4/atomic.m.u; % hopefully = 1
atomic.mass.proton   = 1.007276467/atomic.m.u; 

% hydrogen
atomic.mass.H1   = 1.0078250/atomic.m.u; % abundance 99.9885%
atomic.mass.H2   = 2.0135532/atomic.m.u; % abundance  0.0115%
atomic.mass.H3   = 3.0160492/atomic.m.u; % unstable

% carbon
atomic.mass.C12  =        12/atomic.m.u; % abundance 98.93%
atomic.mass.C13  = 13.003355/atomic.m.u; % abundance  1.109%
atomic.mass.C14  = 14.003241/atomic.m.u; % abundance  1 part per trillion

% nitrogen
atomic.mass.N14  = 14.003074/atomic.m.u; % abundance 99.632%
atomic.mass.N15  = 15.000109/atomic.m.u; % abundance  0.368%

% oxygen
atomic.mass.O16  = 15.994915/atomic.m.u; % abundance 99.762%
atomic.mass.O17  = 16.999132/atomic.m.u; % abundance  0.037%
atomic.mass.O18  = 17.999160/atomic.m.u; % abundance  0.2%

% fluorine
atomic.mass.F19  = 18.998403/atomic.m.u; % abundance  100%

% sodium
atomic.mass.Na23 = 22.989768/atomic.m.u; % abundance  100%

% sulfur
atomic.mass.S32 = 31.9720711/atomic.m.u; % abundance  95%

% chlorine
atomic.mass.Cl35 = 34.968853/atomic.m.u; % abundance 75.78%
atomic.mass.Cl37 = 36.965903/atomic.m.u; % abundance 24.22%

%% Output
% if nargin==0
%     toggle = false;
% end
% if toggle
%     format longE
%     prt.disp ('***************************************************************')
%     prt.disp ('Using atomic units throughout WavePacket!')
%     prt.disp ('***************************************************************')
%     prt.disp ('   ')
%     prt.disp (atomic.m)
%     prt.disp (atomic.l)
%     prt.disp (atomic.E)
%     prt.disp (atomic.T)
%     prt.disp (atomic.w)
%     prt.disp (atomic.t)
%     prt.disp (atomic.o)
%     prt.disp (atomic.d)
%     prt.disp (atomic.F)
%     prt.disp (atomic.I)
%     prt.disp ('***************************************************************')
%     prt.disp ('Most important atomic masses')
%     prt.disp ('***************************************************************')
%     prt.disp ('   ')
%     prt.disp (atomic.mass)
%     format short
% end

