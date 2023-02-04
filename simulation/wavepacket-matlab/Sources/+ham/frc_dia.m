%--------------------------------------------------------------------------
%
% Compute diabatic force matrix at positions r
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function frc_mats = frc_dia(r)
global hamilt space

frc_mats = cell(space.n_dim,1);
for d=1:space.n_dim
    frc_mats{d} = zeros (hamilt.coupling.n_eqs,hamilt.coupling.n_eqs,length(r{1}));
end

for p = 1:hamilt.coupling.n_eqs % diagonal
    for q = p:hamilt.coupling.n_eqs % off-diagonal
        F_pq = F ( hamilt.pot{p,q}, r );
        for d=1:space.n_dim
            frc_mats{d}(p,q,:) = F_pq{d};
            frc_mats{d}(q,p,:) = F_pq{d};
        end
    end
end
end
        
