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

hop_time_counter = zeros(obj.n_p,1);

% Dont allow return hops on the next time step
% Hops on same direction as hop before are still allowed
obj.hop_direction(obj.hop_direction>0) = obj.hop_direction(obj.hop_direction>0) - 1;
obj.hop_direction(obj.hop_direction<0) = obj.hop_direction(obj.hop_direction<0) + 1;

% Reset hop time variable such that a hop in the next time step
% is not falsely not possible
for m = 2:hamilt.coupling.n_eqs
    obj.hop_delta_t{m-1,m}(obj.hop_direction~=0) = obj.hop_delta_t{m-1,m}(obj.hop_direction~=0) + time.steps.s_delta;
    obj.hop_delta_t{m,m-1}(obj.hop_direction~=0) = obj.hop_delta_t{m,m-1}(obj.hop_direction~=0) + time.steps.s_delta;
end

% Current quantum state "m"
for m = 1:hamilt.coupling.n_eqs
    
    % Find indices of trajectories initially in "active" channel "m"
    ind_m = find (obj.cha==m);
    if ~isempty (ind_m)
        
        ind_hop_down=[];
        ind_hop_up=[];
        ind_hop_both_up_first=[];
        ind_hop_both_down_first=[];
        
        if(m~=1)
            % Indices to perform hop to m-1 suggested by LZGD
            ind_hop_down = get_ind_hop(obj, ind_m(obj.hop_direction(ind_m) == 0) , m , m-1);
        end
        
        if(m~=hamilt.coupling.n_eqs)
            % Indices to perform hop to m+1 suggested by LZGD
            ind_hop_up = get_ind_hop(obj, ind_m(obj.hop_direction(ind_m) == 0) , m , m+1);
        end
        
        if ~isempty (ind_hop_up) || ~isempty (ind_hop_down)
            
            % Indicies which cam either hop to m+1 or to m-1
            ind_hop_only_up   = setdiff(ind_hop_up  ,ind_hop_down);
            ind_hop_only_down = setdiff(ind_hop_down,ind_hop_up  );
            
            % Indicies which can hop to m+1 and to m-1
            ind_hop_both = intersect( ind_hop_down , ind_hop_up );
            
            if ~isempty (ind_hop_both)
                % Indicies which came first depending on the
                % step time
                ind_hop_both_up_first   = ind_hop_both( obj.hop_delta_t{m,m+1}(ind_hop_both) >  obj.hop_delta_t{m,m-1}(ind_hop_both) );
                ind_hop_both_down_first = ind_hop_both( obj.hop_delta_t{m,m+1}(ind_hop_both) <= obj.hop_delta_t{m,m-1}(ind_hop_both) );
            end
            
            % Try to perform hop to m+1
            ind_hop_up_first          = union( ind_hop_only_up , ind_hop_both_up_first );
            [ind_hop_up_allowed , ind_hop_up_frustrated] = perform_hop (obj, ind_hop_up_first, m, m+1);
            
            obj.hop_direction(ind_hop_up_allowed) = 2;
            hop_time_counter(ind_hop_up_allowed) = hop_time_counter(ind_hop_up_allowed) + 1;
            ind_hop_both_up_frustrated = intersect( ind_hop_up_frustrated , ind_hop_both ); % frustrated hops can try to hop to lower level
            
            % Try to perform hop to m-1
            ind_hop_down_first_and_up_frustrated = union( union( ind_hop_only_down , ind_hop_both_down_first ) , ind_hop_both_up_frustrated );
            [ind_hop_down_allowed , ind_hop_down_frustrated]    = perform_hop (obj, ind_hop_down_first_and_up_frustrated, m, m-1);
            obj.hop_direction(ind_hop_down_allowed) = -2;
            hop_time_counter(ind_hop_down_allowed) = hop_time_counter(ind_hop_down_allowed) + 1;
            
            % Try to perform hop to m+1 for frustraded hops of
            % hopping down
            ind_hop_both_down_frustated = setdiff( intersect( ind_hop_down_frustrated , ind_hop_both ) , ind_hop_both_up_frustrated);
            ind_hop_up_second_allowed = perform_hop (obj, ind_hop_both_down_frustated, m, m+1);
            obj.hop_direction(ind_hop_up_second_allowed) = 2;
            hop_time_counter(ind_hop_up_second_allowed) = hop_time_counter(ind_hop_up_second_allowed) + 1;
            
        end
        
    end
    
end

% Try to hop further to lower levels
for m = (hamilt.coupling.n_eqs-1):-1:1
    
    % Find indices of trajectories initially in "active" channel "m"
    ind_m = find (obj.cha==m);
    if ~isempty (ind_m)
        
        ind_hop_down=[];
        
        if(m~=1)
            ind_m_direction = ind_m( obj.hop_direction(ind_m) < 0 );
            
            ind_hop_down = get_ind_hop(obj, ind_m_direction , m , m-1);
            
            ind_hop_down = intersect( ind_hop_down , ...
                ind_m( ...
                obj.hop_delta_t{m+1,m}(ind_m) + obj.hop_interval{m+1,m}(ind_m) ...
                >= obj.hop_delta_t{m,m-1}(ind_m) - obj.hop_interval{m,m-1}(ind_m) ...
                ) );
        end
        
        if ~isempty (ind_hop_down)
            % Try to perform hop to m-1
            ind_hop_allowed = perform_hop (obj, ind_hop_down, m, m-1);
            
            hop_time_counter(ind_hop_allowed) = hop_time_counter(ind_hop_allowed) + 1;
        end
    end
end

% Try to hop further to lower levels
for m = 2:hamilt.coupling.n_eqs
    
    % Find indices of trajectories initially in "active" channel "m"
    ind_m = find (obj.cha==m);
    if ~isempty (ind_m)
        
        ind_hop_up=[];
        
        if(m~=hamilt.coupling.n_eqs)
            ind_m_direction = ind_m( obj.hop_direction(ind_m) > 0 );
            
            ind_hop_up = get_ind_hop(obj, ind_m_direction , m , m+1);
            
            ind_hop_up = intersect( ind_hop_up , ...
                ind_m( ...
                obj.hop_delta_t{m-1,m}(ind_m) + obj.hop_interval{m-1,m}(ind_m)...
                >= obj.hop_delta_t{m,m+1}(ind_m) - obj.hop_interval{m,m+1}(ind_m) ...
                ) );
        end
        
        if ~isempty (ind_hop_up)
            % Try to perform hop to m-1
            ind_hop_allowed = perform_hop (obj, ind_hop_up, m, m+1);
            
            hop_time_counter(ind_hop_allowed) = hop_time_counter(ind_hop_allowed) + 1;
        end
    end
end

for m = 1:hamilt.coupling.n_eqs
    time.counter.hops(m) = time.counter.hops(m) + sum( hop_time_counter == m );
end

end