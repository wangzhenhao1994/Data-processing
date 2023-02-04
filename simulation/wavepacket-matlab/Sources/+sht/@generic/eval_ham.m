%--------------------------------------------------------------------------
%
% Evaluate Hamiltonian matrices attached
% to ALL trajectories in ALL channels
% for use in surface hopping methods.
%
% Depending on the surface hopping variant,
% these matrices are used for the following:
% - Calculating energy gaps and/or
% - Integrating the TDSEs attached
%
% In the DIABATIC case, this function provides the following:
%
% trajectories.ham (full potential matrix)
% 
% In the ADIABATIC case, this function provides the following:
%
% trajectories.ham   (diagonal only: adiabatic potential energy surfaces)
% trajectories.nac_1 (off-diagonal only: non-adiabatic coupling, 1st order)
% where the latter ones are optional (for example, not needed in LZ_2)
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu-Araujo
%
% see the README file for license details.

function eval_ham (obj)
global hamilt

if strcmpi(hamilt.coupling.represent,'adi') % Adiabatic representation
    
    % Adiabatic potential energy surfaces
    for m = 1:hamilt.coupling.n_eqs
        obj.ham{m,m}(:) = obj.D_new(m,:);
    end
    
else % Diabatic formulation
    
    for m = 1:hamilt.coupling.n_eqs
        if isa (hamilt.pot{m,m}, 'pot.empty')
            obj.ham{m,m} = zeros(size(obj.pos{1}));
        else
            obj.ham{m,m} = V ( hamilt.pot{m,m}, obj.pos );
        end
        for n = m+1:hamilt.coupling.n_eqs
            if isa (hamilt.pot{m,n}, 'pot.empty')
                V_mn = zeros(size(obj.pos{1}));
            else
                V_mn = V ( hamilt.pot{m,n}, obj.pos );
            end
            obj.ham{m,n} = V_mn;
            obj.ham{n,m} = V_mn;
        end
    end
end

end

