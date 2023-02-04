%--------------------------------------------------------------------------
%
% Update class properties according from level m into level n
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Leonardo Cancissu Araujo
%
% see the README file for license details.

function after_hop (obj,mom_new,ind_h_allowed,allowed,m,n)
global space time

% Inherit from initialization of superclass
after_hop@sht.sssh_gd (obj,mom_new,ind_h_allowed,allowed,m,n);

% Current values
pos = cell(space.n_dim,1);
mom = cell(space.n_dim,1);
frc = cell(space.n_dim,1);
for d=1:space.n_dim
    pos{d} = obj.pos{d}(ind_h_allowed);
    mom{d} = obj.mom{d}(ind_h_allowed);
    frc{d} = obj.frc{d}(ind_h_allowed);
end

% Perform verlet two times backwards
[pos_old, mom_old, frc_old, pot_mats_old, frc_mats_old, D_old, U_old] ...
    = verlet_prop_back(obj,pos,mom,frc,n);

[pos_oldold, mom_oldold, frc_oldold, pot_mats_oldold, frc_mats_oldold, D_oldold, U_oldold] ...
    = verlet_prop_back(obj,pos_old,mom_old,frc_old,n);

% Update previous variables
obj.pot_mat_old   (:,:,ind_h_allowed) = pot_mats_old;
obj.pot_mat_oldold(:,:,ind_h_allowed) = pot_mats_oldold;

obj.U_old   (:,:,ind_h_allowed) = U_old;
obj.U_oldold(:,:,ind_h_allowed) = U_oldold;

obj.D_old   (:,ind_h_allowed) = D_old;
obj.D_oldold(:,ind_h_allowed) = D_oldold;

for d=1:space.n_dim
    obj.pos_old   {d}(ind_h_allowed) = pos_old   {d};
    obj.pos_oldold{d}(ind_h_allowed) = pos_oldold{d};
    
    obj.mom_old   {d}(ind_h_allowed) = mom_old   {d};
    obj.mom_oldold{d}(ind_h_allowed) = mom_oldold{d};
    
    obj.frc_old   {d}(ind_h_allowed) = frc_old   {d};
    obj.frc_oldold{d}(ind_h_allowed) = frc_oldold{d};
    
    obj.frc_mat_old   {d}(:,:,ind_h_allowed) = frc_mats_old   {d};
    obj.frc_mat_oldold{d}(:,:,ind_h_allowed) = frc_mats_oldold{d};
end

% Counting the extra diagonalizations produced by the gradient descent method
time.counter.diag.correction = time.counter.diag.correction + 2 * length(ind_h_allowed);

end

