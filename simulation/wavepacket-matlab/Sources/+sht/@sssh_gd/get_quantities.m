%--------------------------------------------------------------------------
%
% Give quantities needed to check momentum rescaling and to perform hop
% The variables are of the size of ind_hop
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function [pos,mom,pot_mats,frc_mats,U,D] = get_quantities(obj,ind,m,n)
global space

pos      = cell(space.n_dim,1);
mom      = cell(space.n_dim,1);
frc_mats = cell(space.n_dim,1);

for d = 1:space.n_dim
    % Position/momentum while hopping
    pos{d}      = obj.hop_pos    {m,n}{d}(ind);
    mom{d}      = obj.hop_mom    {m,n}{d}(ind);
    frc_mats{d} = obj.hop_frc_mat{m,n}{d}(:,:,ind);
end

pot_mats = obj.hop_pot_mat{m,n}(:,:,ind);
U        = obj.hop_U      {m,n}(:,:,ind);
D        = obj.hop_D      {m,n}(:,  ind);
end

