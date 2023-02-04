%--------------------------------------------------------------------------
%
% Perform the trajectory hopping
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function traj_hop (obj)
global hamilt time

ind_cha = cell(hamilt.coupling.n_eqs,1);

% Save current indicies of trajectories which are on the same level
for m = 1:hamilt.coupling.n_eqs
    ind_cha{m} = find (obj.cha==m);
end

% Current quantum state "m"
for m = 1:hamilt.coupling.n_eqs
    
    % Find indices of trajectories initially in "active" channel "m"
    ind_m = ind_cha{m};
    if ~isempty (ind_m)
        
        % Uniform random numbers in the interval (0,1)
        zeta = rand(size(ind_m));
        
        % Initialize summations of hopping probabilities
        summ = zeros(size(ind_m));
        
        % Loop over all other quantum states "n"
        for n = 1:hamilt.coupling.n_eqs
            
            % Adiabatic representation: neighboring states only
            % Diabatic representation: All other states
            if n~=m && (n==m+1 || n==m-1 || strcmpi(hamilt.coupling.represent,'dia') || ~obj.consec_hops)
                
                % Get probabilities of hopping (from sub-classes)
                probable = prob_hop (obj,m,n,ind_m);
                
                % If probabilities are negative, set them to zero
                probable ( probable<0 ) = 0;
                
                % Summing up probabilities, see e.g. Eq. (10) of doi:10.1063/1.5000843
                prev = summ;
                summ = summ + probable;
                
                % Find indices of hopping trajectories by comparison with
                % zeta, a uniform random number in the interval (0,1)
                ind_hop = ind_m (zeta>prev & zeta<summ );
                
                % Check momentum rescaling and perform hop
                if ~isempty(ind_hop)
                    
                    % Get quantities needed to check momentum rescaling and to perform hop
                    % The variables are of the size of ind_hop
                    [pos,mom,pot_mats,frc_mats,U,D] = get_quantities(obj,ind_hop,m,n);
                    
                    % Energy conservation; allowed and frustated hops
                    if obj.rescale
                        % Check momentum rescaling for proposed trajectories
                        [mom_new, allow, forbid] = mom_rescaling(obj,mom,pot_mats,frc_mats,U,D,m,n);
                    else
                        % If no momentum rescaling, then let mom unchanged and allow all proposed hops
                        mom_new   = mom;
                        allow  = true  ( length(ind_hop) , 1 );
                        forbid = false ( length(ind_hop) , 1 );
                    end
                    
                    % Save indicies of trajectories which are allowed / forbidden
                    ind_hop_allowed   = ind_hop(allow);
                    ind_hop_forbidden = ind_hop(forbid);
                    
                    % Perform hop for allowed trajectories
                    if ~isempty(ind_hop_allowed)
                        % Update class properties according to the new level n
                        after_hop (obj,mom_new,ind_hop_allowed,allow,m,n)
                    end
                    
                    % Store ind_hop_allowed and ind_hop_forbidden hops in global variables
                    % Needed for scatter plots
                    store_allowed_forbidden (obj, ind_hop_allowed, ind_hop_forbidden, m, n)
                    
                    % Display allowed and forbidden hops
                    % including relevant quantities
                    if obj.verbose
                        disp_hop(obj,pos,mom,mom_new,D,ind_hop_allowed,ind_hop_forbidden,allow,m,n)
                    end
                    
                    time.counter.hops(abs(m-n)) = time.counter.hops(abs(m-n)) + length(ind_hop_allowed);
                    
                end
                
            end
        end
    end
end
end