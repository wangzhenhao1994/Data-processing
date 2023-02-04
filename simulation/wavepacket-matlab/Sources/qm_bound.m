%------------------------------------------------------------------------------
%
% Solves the time-independent Schroedinger equation to 
% obtain the bound states of a suitable quantum system.
%
%------------------------------------------------------------------------------


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-20.. Burkhard Schmidt
%               2007-2011 Ulf Lorenz
%
% see the README file for license details.

function qm_bound ()

% Main variables are global throughout;
global hamilt state

% Initializes general information and sets up log files.
prt.init (mfilename('fullpath'));

% Bound state calculations for wavefunctions only
if ~isa(state,'wave') 
    prt.error ('Bound state calculations for wavefunctions only')
end

% Initialize spatial discretization for each degree of freedom
dof.init (state);

% Initialize Hamiltonian operator
init_ham (state);

% Solve eigenproblem for Hamiltonian matrix
init  ( hamilt.eigen );
disp  ( hamilt.eigen );
setup ( hamilt.eigen );
symm  ( hamilt.eigen );
diag  ( hamilt.eigen );

% Initialize expectation values and uncertainties of observables
obs.init;

%% Beginning of main loop over bound states
for step = 1 : hamilt.eigen.number
    
    % Extract psi from columns(!) of eigenvector matrix and normalize
    eigen ( state, step );
    
    % Transform to adiabatic representation (if desired)
    adiabatic ( state, step, 'dia2adi' ); % needs more testing ?!?
    
    % Expectation values and uncertainties of observables
    observe ( state, step );
    
    % Logging and plot title
    obs.log ( step );       
             
    % Show visualization of densities and expectation values
    vis.show ( state, step );

    % save the wave function
    save ( state,  step );
    
% End of main loop     
end

% Output clock/date/time
prt.clock;

end
