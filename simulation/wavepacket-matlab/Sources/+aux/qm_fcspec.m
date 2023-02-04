%------------------------------------------------------------------------------
%
% Franck-Condon spectrum
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-20xy Burkhard Schmidt
%
% see the README file for license details.

function qm_fcspec ( state )

global expect hamilt;

% Initializes general information and sets up log files.
prt.init (mfilename('fullpath'));

% Initialize spatial discretization
dof.init (state);

% Initialize Hamiltonian operator
init_ham (state);

% Initialize wave functions / densities
init_obj (state);

% Initialize expectation values and uncertainties of observables
obs.init;

% Expectation values and uncertainties
observe ( state, 1 );

% Load the general data of the saved calculation.
% Sets most global variables, except for: plots
% which should be set before running this function
load_0 (state, true);

%% Main loop over time steps 
intensity = zeros(hamilt.eigen.stop-hamilt.eigen.start+1,1);
for step = hamilt.eigen.start:hamilt.eigen.stop
    
    % Load states from data files
    load ( state, step+1 );

    % Expectation values and uncertainties
    observe ( state, step+1 );
    
    % Spectral intensity as the square of the overlap
    intensity(step+1) = abs(wave.braket(state.ini,state.dvr))^2;
    
end

% Output clock/date/time
prt.clock;

figure(7);clf
stem(expect.total',intensity)
xlabel('energy')
ylabel('intensity')

end
