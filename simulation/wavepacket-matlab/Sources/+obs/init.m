%------------------------------------------------------------------------------
%
% Initialize expectation values and uncertainties of basic observables
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% see the README file for license details.

function init
global expect hamilt space state time

switch(lower(class(state)))
    
    case 'ket'
        
        state.y = zeros ( time.steps.m_number, length(state.D));
        
    case 'rho'
        
        state.y = zeros ( time.steps.m_number, length(state.C));
        
    otherwise
        
        % Population threshold for logging, plotting
        expect.min_pop = 10^-3;
        
        % Population
        expect.pop = obs.generic ('pop');
        
        % Additional multiplicative operators
        if isfield(hamilt, 'amo')
            for p = 1:length(hamilt.amo)
                if ~isempty (hamilt.amo{p})
                    expect.amo{p} = obs.generic('amo');
                    expect.amo{p}.ind = p;
                end
            end
        end
        
        % DVR and FBR for each spatial dimension
        for k = 1:space.n_dim
            expect.pos{k} = obs.generic('pos');
            expect.pos{k}.ind = k;
            expect.mom{k} = obs.generic('mom');
            expect.mom{k}.ind = k;
        end
        
        % Potential/kinetic energy
        expect.pot = obs.generic ('pot');
        expect.kin = obs.generic ('kin');
        
        % Total energy
        expect.total = zeros ( time.steps.m_number, 1 );
        
end
