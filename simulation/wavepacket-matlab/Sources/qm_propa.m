%------------------------------------------------------------------------------
%
% If "state" is an object of class "wave":
% Fully quantum mechanical dynamics:
% Solves the time-dependent Schroedinger equation 
% by propagating (coupled) wave function(s) in time.
%
% If "state" is an object of class "ket":
% Fully quantum mechanical dynamics:
% Solves the time-dependent Schroedinger equation
% by propagating state vectors in time.
%
% If "state" is an object of class "rho":
% Fully quantum mechanical dynamics:
% Solves the quantum Liouville-von Neumann equation
% by propagating state vectors in time.
%
% If "state" is an object of class "traj":
% Fully classical dynamics:
% Solves the classical Liouville equation by 
% propagating phase-space density(s) in time.
%
% If "state" is an object of class "mssh", "fssh", or "sssh":
% Mixed quantum-mechanical / classical dynamics:
% Solves the quantum-classical Liouville equation by 
% propagating (coupled) phase-space density(s) in time.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-20.. Burkhard Schmidt
%               2007-2011 Ulf Lorenz
%
% see the README file for license details.

function qm_propa (string1, param1, param2)

%% Choice of propagators
global state time

if nargin<2
    param1=[];
end

if nargin<3
    param2=[];
end

if isa(state,'wave')                     % quantum wavefunctions, represented on grids
    if nargin==0
        string1 = 'strang';              % Default for choice of propagator
    end
    switch lower(string1)
        case 'diff_2'                    % Second order differencing
            time.propa = tmp.wave.diff_2;
        case {'lie','trotter'}           % First order Lie-Trotter splitting
            time.propa = tmp.wave.trotter;
        case {'strang', 'marchuk'}       % Second order Strang-Marchuk splitting
            time.propa = tmp.wave.strang;
        case 'cheby_real'                % Chebychev polynomials in real time
            time.propa = tmp.wave.cheby_real (param1, param2);
        case 'cheby_imag'                % Chebychev polynomials in imaginary time
            time.propa = tmp.wave.cheby_imag (param1, param2);
        otherwise
            prt.error ('Wrong choice of propagator for wavefunctions')
    end
    
elseif isa(state,'traj')                 % classical and (hopping!) trajectories
    if nargin==0
        string1 = 'leap_frog';           % Default for choice of propagator
    end
    switch lower(string1)
        case 's_verlet'                  % Stoermer-Verlet
            time.propa = tmp.traj.propagators.s_verlet;
        case 'v_verlet'                  % Velocity Verlet
            time.propa = tmp.traj.propagators.v_verlet;
        case {'march_frog','lie','trotter'}    % March frog (Lie-Trotter)
            time.propa = tmp.traj.propagators.march_frog;
        case {'leap_frog','strang','marchuk'}  % Leap frog (Strang-Marchuk)
            time.propa = tmp.traj.propagators.leap_frog;
        case 'beeman'                    % Beeman 3rd order
            time.propa = tmp.traj.propagators.beeman;
        case 'yoshida'                   % Yoshida 4th order
            time.propa = tmp.traj.propagators.yoshida;
        otherwise
            prt.error ('Wrong choice of propagator for trajectories')
    end
    
elseif isa(state,'ket')  % this includes sub-class 'rho'
    if nargin==0
        string1 = 'matlab';              % Matlab built-in ODE integrators
    end
    switch lower(string1)
        case 'matlab'                    % Matlab built-in ODE integrators
            time.propa = tmp.ode.matlab (param1, param2);
        otherwise
            prt.error ('Wrong choice of propagator for ket//rho objects')
    end
    
else
    prt.error ('Wrong choice of WavePacket main class')
end

%% Initialization

% Initializes general information and sets up log files.
prt.init (mfilename('fullpath'));

% Initialize spatial discretization for each degree of freedom
dof.init (state);

% Initialize Hamiltonian operator
init_ham (state);

% Initialize temporal discretization
init (time.steps);

% Initialize the electric field
efi.init;

% Initialize wave functions / densities
init_obj (state);

% Transform to adiabatic representation (if desired)
adiabatic ( state, -1 );

% Initialize expectation values and uncertainties of observables
obs.init;

%% Main loop over time steps (step=1 is the initial step)
for step = 1 : time.steps.m_number
    
    % Numerical propagation using pde solvers, possibly with absorbing boundary conditions
    propagate ( state, step );
    
    % Transform to adiabatic representation (if desired)
    adiabatic ( state, step, 'dia2adi' );
    
    % Expectation values and uncertainties of observables
    observe ( state, step );
    
    % Logging and plot title
    obs.log ( step );
        
    % Show visualization of densities and expectation values
    vis.show ( state, step );
    
    % Transform back to diabatic representation
    adiabatic ( state, step, 'adi2dia' );
        
    % Save the states in data files
    save ( state, step );
    
    % End of main loop
end

% Output clock/date/time
prt.clock;

end
