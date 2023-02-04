%--------------------------------------------------------------------------
%
% Evaluate potential energies and forces which
% are needed for all types of fully classical
% or quantum-classical trajectory propagations.
% 
% Using either diabatic or adiabatic representation.
% In the latter case, all the diabatic potential and
% force matrices are not only used for calculating 
% adiabatic potentials and forces but they are also
% saved for use in tdse_ham.m (SH and FSSH only)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu-Araujo
%
% see the README file for license details.

function eval_V_F (obj, state, step_size, init)
global hamilt space time

% Adiabatic representation
if(strcmpi(hamilt.coupling.represent,'adi'))
    
    % Calculate and save all diabatic potentials in advance
    state.pot_mat = ham.pot_dia(state.pos);

    % Compute adiabatic potential matrix and eigenvector matrix
    U_old = state.U_new;
    [state.U_new , state.D_new] = ham.pot_eig_adi(state.pot_mat);

    % Ensure continuity of eigenvector matrix
    if( hamilt.coupling.n_eqs > 2 )
        ensure_continuity_eig (state, U_old);
    end

    % Calculate and save all diabatic forces in advance
    state.frc_mat = ham.frc_dia(state.pos);

    if(~init)
        solve_tdse (obj, state, step_size, state.D_new, state.U_new, U_old) % adiabatic
    end

    % Loop over (coupled) channels m
    for m = 1:hamilt.coupling.n_eqs

        % Indices of trajectories in that channel
        ind_m = find(state.cha==m);

        if(~isempty(ind_m))

            % Save adiabatic potential energies
            state.pot(ind_m) = state.D_new(m,ind_m);

            frc_mats_ind_m = cell(space.n_dim,1);

            for d=1:space.n_dim
                frc_mats_ind_m{d} = state.frc_mat{d}(:,:,ind_m);
            end

            psi_ind_m = cell(hamilt.coupling.n_eqs,1);

            for n=1:hamilt.coupling.n_eqs
                psi_ind_m{n} = state.psi{n} (ind_m);
            end

            frc_adi = eq_frc (time.eq_motion, state.pot_mat(:,:,ind_m) , frc_mats_ind_m,...
                state.D_new(:,ind_m), state.U_new(:,:,ind_m), m, psi_ind_m);

            for d=1:space.n_dim
                state.frc{d}(ind_m) = frc_adi{d};
            end
        end
    end

    time.counter.diag.dynamics = time.counter.diag.dynamics + state.n_p;
    
% Diabatic representation
else
    
    % Loop over coupled channels
    r_m = cell(space.n_dim,1);
    for m = 1:hamilt.coupling.n_eqs
        
        % Get position vectors of trajectories in the respective channel
        for d=1:space.n_dim
            r_m{d} = state.pos{d}(state.cha==m);
        end
        
        % If there are any trajectories in the respective channel
        if ~isempty(r_m)
            
            % Get corresponding potentials and forces
            if isa ( hamilt.pot{m,m}, 'pot.empty' ) % free particle
                V_m = zeros(size(r_m{1}));
                F_m = cell(space.n_dim,1);
                for d=1:space.n_dim
                    F_m{d} = zeros(size(r_m{1}));
                end
            else
                V_m = V ( hamilt.pot{m,m}, r_m );
                F_m = F ( hamilt.pot{m,m}, r_m );
            end
            
            % Save corresponding potentials and forces
            state.pot (state.cha==m) = V_m;
            for d=1:space.n_dim
                state.frc{d}(state.cha==m) = F_m{d};
            end
            
            if (state.extra_tdse_solving)
                tdse_dia (state) % diabatic
            end
        end
        
    end
    
end

end