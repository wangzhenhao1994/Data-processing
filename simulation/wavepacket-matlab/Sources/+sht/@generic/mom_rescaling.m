%--------------------------------------------------------------------------
%
% Rescaling momenta of classical particles to enforce
% energy conservation in surface hopping trajectories
%
% Note that upon hopping from lower to upper states
% it may be that the momenta are not sufficient in
% which case the energy cannot be conserved, hence
% the hopping is not performed ("frustrated hops")
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function [mom_new, allowed, forbidden] = mom_rescaling(obj,mom,pot_mats,frc_mats,U,D,m,n)

% Choose vector for momentum rescaling
% Maybe sca_nac should be renamed
switch obj.sca_nac
    
    % Along momenta
    case 0
        vec = mom;
        
        % Along NAC vector
    case 1
        vec = ham.nac_mn(pot_mats,frc_mats,U,D,m,n);
end

% Perform momentum rescaling along chosen vector
[mom_new, allowed, forbidden] = mom_rescaling_vec (obj,mom,D,vec,m,n);

end
 