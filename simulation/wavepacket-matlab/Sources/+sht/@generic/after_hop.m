%--------------------------------------------------------------------------
%
% Update class properties according to the new level n
% 
% One dummy arguments are to ensure consistency with 
% equally named method in SSSH_GD class definition
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function after_hop (obj,mom_new,ind_h_allowed,allowed,~,n)
global space hamilt

% Perform hopping
obj.cha(ind_h_allowed) = n;

% Momentum rescaling
for d=1:space.n_dim
    obj.mom{d}(ind_h_allowed) = mom_new{d}(allowed);
end

if strcmpi(hamilt.coupling.represent,'adi') % Adiabatic representation
    % New potential energy
    obj.pot(ind_h_allowed) = obj.D_new(n,allowed)';
    
    % Calculation of new force:
    frc_mats_allowed = cell(space.n_dim,1);
    
    for d=1:space.n_dim
        frc_mats_allowed{d} = obj.frc_mat{d}(:,:,allowed);
    end
    
    frc_adi = ham.frc_adi ( obj.pot_mat(:,:,allowed) , ...
        frc_mats_allowed , obj.U_new(:,:,allowed) , n );
    
    % Store new force
    for d=1:space.n_dim
        obj.frc{d}(ind_h_allowed) = frc_adi{d};
    end

else % Diabatic representation
    % New potential energy
    obj.pot(ind_h_allowed) = obj.ham{n,n}(allowed);
    
    % FORCE NEEDS TO BE UPDATED...
end

end

