%--------------------------------------------------------------------------
%
% Update class properties according from level m into level n
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function after_hop (obj,mom_new,ind_h_allowed,allowed,m,n)
global space time

% Perform hopping
obj.cha(ind_h_allowed) = n;

% Save pos, mom and frc of last step of GD as current time step
% -> mom is the new rescaled momentum
for d=1:space.n_dim
    obj.pos{d}(ind_h_allowed) = obj.hop_pos  {m,n}{d}(ind_h_allowed);
    obj.frc{d}(ind_h_allowed) = obj.hop_frc_n{m,n}{d}(ind_h_allowed);
    
    obj.mom{d}(ind_h_allowed) = mom_new{d}(allowed);
end

% Integrate to true current time step (Verlet):

% Propagate momenta by half time step
for d = 1:space.n_dim
    obj.mom{d}(ind_h_allowed) = obj.mom{d}(ind_h_allowed) + ...
        0.5 * obj.hop_delta_t{m,n}(ind_h_allowed) .* obj.frc{d}(ind_h_allowed);
end

% Propagate positions by full time step
for d = 1:space.n_dim
    obj.pos{d}(ind_h_allowed) = obj.pos{d}(ind_h_allowed) + ...
        obj.hop_delta_t{m,n}(ind_h_allowed) .* obj.mom{d}(ind_h_allowed) / space.dof{d}.mass;
end

% Computing adiabatic force:
pos_allowed = cell(space.n_dim,1);
for d=1:space.n_dim
    pos_allowed{d} = obj.pos{d}(ind_h_allowed);
end

% Calculate and save all diabatic potentials in advance
pot_mats_allowed = ham.pot_dia(pos_allowed);

% Compute adiabatic potential matrix and eigenvector matrix
[U_1 , D_1] = ham.pot_eig_adi(pot_mats_allowed);
obj.D_new(:,ind_h_allowed)   = D_1;
obj.U_new(:,:,ind_h_allowed) = U_1;

% New potential energy
obj.pot(ind_h_allowed) = D_1(n,:);

% Calculate and save all diabatic forces in advance
frc_mats_allowed = ham.frc_dia(pos_allowed);

frc_n = ham.frc_adi ( pot_mats_allowed, frc_mats_allowed , U_1 , n );
for d=1:space.n_dim
    obj.frc{d}(ind_h_allowed) = frc_n{d};
end

% Propagate momenta by half time step
for d = 1:space.n_dim
    obj.mom{d}(ind_h_allowed) = obj.mom{d}(ind_h_allowed) + ...
        0.5 * obj.hop_delta_t{m,n}(ind_h_allowed) .* obj.frc{d}(ind_h_allowed);
end

% Counting the extra diagonalizations produced by the gradient descent method
time.counter.diag.correction = time.counter.diag.correction + length(ind_h_allowed);

%update variables
obj.pot_mat(:,:,ind_h_allowed) = pot_mats_allowed;
for d = 1:space.n_dim
    obj.frc_mat{d}(:,:,ind_h_allowed) = frc_mats_allowed{d};
end

end
