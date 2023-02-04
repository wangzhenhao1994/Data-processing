%------------------------------------------------------------------------------
%
% Regenerates some global data needed to restart a saved calculation.
%
% The problem: If you save a calculation and this calculation is interrupted,
% e.g., by a crash, then all the global variables are not written out. They are
% needed, however, to load the wave functions again (wave.load() is written in
% such a way). This function regenerates the global variables, so that the
% calculation can later be started by the rerun() function.
%
% To use this function, setup the calculation, but then instead of qm_propa() et al.
% call regenerate().
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt
%               2007-2009 Ulf Lorenz
%               2008 Burkhard Schmidt
%               2011 Ulf Lorenz
%
% see the README file for license details.

function qm_regenerate()

global time;

% Initializes general information and sets up log files.
prt.init (mfilename('fullpath'));

%% Initialize coupling scheme and spatial/temporal discretization
init.grid;

% Initialize the electric field
init.efield;

% Initialize Hamiltonian operator 
init.hamilt;

% Initialize temporal discretization
tmp.init;

% Initialize wave function
wave.init;

% Initialize expectation values and uncertainties of observables
obs.init;

%% Beginning of main loop over time steps (step=1 is the initial step )
step = 1;
    
    % Numerical propagation using pde solvers, possibly with absorbing boundary conditions
    wave.propa ( step );

    % Transform to adiabatic representation (if desired)
    wave.adiabatic ( 'dia2adi' );
    
    % Expectation values and uncertainties
    wave.observe ( step );
    
    % Logging and plot title
    obs.log ( step );       

    % Get spectrum as Fourier transform of autocorrelation
    wave.spectrum ( step );

    % Plot densities and expectation values
    plot.ket ( step );

    % Transform back to diabatic representation
    wave.adiabatic ( 'adi2dia' );

    % Store the wave function.
    wave.save ( step );
    

% This should save the missing global variables
wave.save( time.steps.m_number );


% Output clock/date/time
prt.clock;
