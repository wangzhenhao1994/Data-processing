%--------------------------------------------------------------------------
%
% Display details of surface hopping
% Function not in used anymore and is 
% probably not working correctly
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function disp_hop(~,pos_cell,mom_cell,mom_new_cell,D,ind_h_allowed,ind_h_forbidden,allowed,m,n)
global space

n_hop   = length(pos_cell{1});

pos     = zeros(n_hop , space.n_dim);
mom     = zeros(n_hop , space.n_dim);
mom_new = nan  (n_hop , space.n_dim);

for d = 1:space.n_dim
    % Position/momentum while hopping
    pos(:,d) = pos_cell{d};
    mom(:,d) = mom_cell{d};
    
    mom_new(allowed,d) = mom_new_cell{d}(allowed);
end

% Potential energy before|after hopping
pot_m = D(m,:)';
pot_n = D(n,:)';

% Kinetic energy before|after hopping from energy conservation
kin_m = zeros(n_hop,1);
for d = 1:space.n_dim
    kin_m = kin_m + ( mom(:,d).^2 ./ (2*space.dof{d}.mass));
end

total = pot_m + kin_m;

kin_n = total - pot_n;

prt.disp (['Hopping from state ' int2str(m) ...
    ' to ' int2str(n) ...
    ': ' int2str(n_hop) ...
    ' trajectories'])

j_allowed = 0;

for j=1:n_hop
    if allowed(j)
        j_allowed = j_allowed + 1;
        doit = 'allowed   ';
        ind_j = ind_h_allowed (j_allowed);
    else
        doit = 'FORBIDDEN ';
        ind_j = ind_h_forbidden (j - j_allowed);
    end
    prt.disp ( [ doit...
        'traj_ind: ' int2str(ind_j)             ': '  ...
        'pos: '      num2str(pos (j,:))         ' | ' ...
        'mom_prev: ' num2str(mom (j,:))         ' | ' ...
        'mom_new: '  num2str(mom_new (j,:))     ' # ' ...
        'pot_m: '    num2str(pot_m (j))         ' | ' ...
        'pot_n: '    num2str(pot_n (j))         ' | ' ...
        'kin_m: '    num2str(kin_m (j))         ' | ' ...
        'kin_n: '    num2str(kin_n (j))         ' | ' ...
        'E_tot: '    num2str(total (j)) ] )
end

prt.disp ('   ')
end
