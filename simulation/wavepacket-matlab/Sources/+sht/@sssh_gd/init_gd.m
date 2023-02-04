%--------------------------------------------------------------------------
%
% Initialization of variables needed for the gradient descent SSSH
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-... Leonardo Cancissu Araujo
%
% see the README file for license details.

function init_gd(obj)
global hamilt space

if(isempty(obj.max_round_gd))
    obj.max_round_gd = 20;
end
if(isempty(obj.acc_gd))
    obj.acc_gd = 10^(-4);
end

obj.hop_delta_t = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs);
obj.hop_interval= cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs);
obj.hop_pos     = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs);
obj.hop_mom     = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs);
obj.hop_ham     = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs);
obj.hop_frc_n   = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs);
obj.hop_pot_mat = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs);
obj.hop_U       = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs);
obj.hop_D       = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs);
obj.hop_frc_mat = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs);

for m=1:hamilt.coupling.n_eqs
    for n=1:hamilt.coupling.n_eqs
        obj.hop_delta_t{m,n} = zeros(obj.n_p,1);
        obj.hop_interval{m,n}= zeros(obj.n_p,1);
        obj.hop_pot_mat{m,n} = zeros (hamilt.coupling.n_eqs,hamilt.coupling.n_eqs,obj.n_p);
        obj.hop_U{m,n}       = zeros (hamilt.coupling.n_eqs,hamilt.coupling.n_eqs,obj.n_p);
        obj.hop_D{m,n}       = zeros (hamilt.coupling.n_eqs,obj.n_p);
        
        obj.hop_pos{m,n}     = cell(space.n_dim,1);
        obj.hop_mom{m,n}     = cell(space.n_dim,1);
        obj.hop_frc_n{m,n}   = cell(space.n_dim,1);
        obj.hop_frc_mat{m,n} = cell(space.n_dim,1);
        
        for d = 1:space.n_dim
            obj.hop_pos{m,n}{d}   = zeros(obj.n_p,1);
            obj.hop_mom{m,n}{d}   = zeros(obj.n_p,1);
            obj.hop_frc_n{m,n}{d} = zeros(obj.n_p,1);
            
            obj.hop_frc_mat{m,n}{d} = zeros (hamilt.coupling.n_eqs,hamilt.coupling.n_eqs,obj.n_p);
        end
        
        obj.hop_ham{m,n} = cell(hamilt.coupling.n_eqs,hamilt.coupling.n_eqs);
        
        for p=1:hamilt.coupling.n_eqs
            for q=1:hamilt.coupling.n_eqs
                obj.hop_ham{m,n}{p,q} = zeros(obj.n_p,1);
            end
        end
    end
end

end
