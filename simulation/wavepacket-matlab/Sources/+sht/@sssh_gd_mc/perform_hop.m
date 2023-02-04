%--------------------------------------------------------------------------
%
% Perform hopping from old level m into new level n
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-.... Leonardo Cancissu Araujo
%
% see the README file for license details.

function [ind_hop_allowed,ind_hop_forbidden] = perform_hop (obj, ind_hop, m, n)
ind_hop_allowed  =[];
ind_hop_forbidden=[];

% Check momentum rescaling and perform hop
if ~isempty(ind_hop)
    
    % Get quantities needed to check momentum rescaling and to perform hop
    % The variables are of the size of ind_hop
    [pos,mom,pot_mats,frc_mats,U,D] = get_quantities (obj,ind_hop,m,n);
    
    % Energy conservation; allowed and frustated hops
    if obj.rescale
        % Check momentum rescaling for proposed trajectories
        [mom_new, allowed, forbidden] = mom_rescaling(obj,mom,pot_mats,frc_mats,U,D,m,n);
    else
        % If no momentum rescaling, then let mom unchanged and allow all proposed hops
        mom_new   = mom;
        allowed   = true  ( length(ind_hop) , 1 );
        forbidden = false ( length(ind_hop) , 1 );
    end
    
    % Save indicies of trajectories which are allowed / forbidden
    ind_hop_allowed   = ind_hop(allowed);
    ind_hop_forbidden = ind_hop(forbidden);
    
    % Perform hop for allowed trajectories
    if ~isempty(ind_hop_allowed)
        % Update class properties according to the new level n
        after_hop (obj,mom_new,ind_hop_allowed,allowed,m,n)
    end
    
    % Store ind_hop_allowed and ind_hop_forbidden hops in global variables
    % Needed for scatter plots
    store_allowed_forbidden (obj, ind_hop_allowed, ind_hop_forbidden, m, n)
    
    % Display allowed and forbidden hops
    % including relevant quantities
    if obj.verbose
        disp_hop(obj.pos,mom,mom_new,D,ind_hop_allowed,ind_hop_forbidden,allowed,m,n)
    end
end
end
