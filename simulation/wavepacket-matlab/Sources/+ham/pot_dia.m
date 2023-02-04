%--------------------------------------------------------------------------
%
% Compute diabatic potential matrix at positions r
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function pot_mats = pot_dia(r)
global hamilt

pot_mats = zeros ( hamilt.coupling.n_eqs , hamilt.coupling.n_eqs , length(r{1}) );

for p = 1:hamilt.coupling.n_eqs % diagonal
    for q = p:hamilt.coupling.n_eqs % off-diagonal
        if isa (hamilt.pot{p,q}, 'pot.empty')
            V_pq = zeros(size(r{1}));
        else
            V_pq = V ( hamilt.pot{p,q}, r );
        end
        pot_mats(p,q,:) = V_pq;
        pot_mats(q,p,:) = V_pq;
    end
end
end
