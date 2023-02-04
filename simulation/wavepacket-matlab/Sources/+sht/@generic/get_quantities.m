%--------------------------------------------------------------------------
%
% Get quantities needed to check momentum rescaling and to perform hop
% The variables are of the size of input argument "ind"
% 
% Two dummy arguments are to ensure consistency with 
% equally named method in SSSH_GD class definition
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function [pos,mom,pot_mats,frc_mats,U,D] = get_quantities(obj,ind,~,~)
global space hamilt

pos      = cell(space.n_dim,1);
mom      = cell(space.n_dim,1);

for d = 1:space.n_dim
    % Position/momentum while hopping
    pos{d}      = obj.pos{d}(ind);
    mom{d}      = obj.mom{d}(ind);
end

% Adiabatic representation
if(strcmpi(hamilt.coupling.represent,'adi'))

    pot_mats = obj.pot_mat(:,:,ind);
    U        = obj.U_new  (:,:,ind);
    D        = obj.D_new  (:,  ind);

    frc_mats = cell(space.n_dim,1);
    for d = 1:space.n_dim
        % Diabatic forces while hopping
        frc_mats{d} = obj.frc_mat{d}(:,:,ind);
    end

% Diabatic representation
else
    D = zeros  (hamilt.coupling.n_eqs,  length(ind));
    
    for m = 1:hamilt.coupling.n_eqs
        D(m,:) = obj.ham{m,m} (ind);
    end
    
    pot_mats = [];
    U        = [];
    frc_mats = [];
end

end
