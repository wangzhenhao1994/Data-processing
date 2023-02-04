%--------------------------------------------------------------------------
%
% Quantum trajectories in adiabatic representation
% with non-classical contributions to forces
%
% see work by Craig C. Martens
% DOI:10.1021/acs.jpca.8b10487
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2020-.... Leonardo Cancissu Araujo
%
% see the README file for license details.
%  comment added by CCM on May 20, 2020

classdef quantraj_adi < tmp.traj.eq_motion.classical & handle
    
    properties (Access = public)
    end    

    methods (Access = public)
            
        % Initialize
        function eq_init (obj, state, propa_type) 
            global space
            
            eq_init@tmp.traj.eq_motion.classical ( obj, state, propa_type)
            
            state.mom_kin = cell(space.n_dim,1);
            
            mom_diff = get_mom_diff (obj,state.pot_mat,state.frc_mat,state.D_new,state.U_new,state.psi);
            
            % Get initial kinetic momentum from "normal" momentum
            for d = 1:space.n_dim
                state.mom_kin{d} = state.mom{d} - mom_diff{d};
            end
        end
        
        % Propagate
        function eq_propa (obj, state, propa_type) 
            global space
            
            state.mom = state.mom_kin;
            
            eq_propa@tmp.traj.eq_motion.classical ( obj, state, propa_type)
                        
            mom_diff = get_mom_diff (obj,state.pot_mat,state.frc_mat,state.D_new,state.U_new,state.psi);
            
            state.mom_kin = state.mom;
                  
            % Get "normal" momentum from kinetic momentum
            for d = 1:space.n_dim
                state.mom{d} = state.mom_kin{d} + mom_diff{d};
            end
            
        end
        
        function solve_tdse(obj, state, step_size, D, U, U_old)
            % Always solve TDSE for quantraj
            tdse_adi(state, step_size, D, U, U_old)
        end
        
        % Propagate
        function frc_adi = eq_frc (obj, pot_mats, frc_mats, D, U, m, psi) 
            
            % Get classical adiabatic force
            frc_adi = eq_frc@tmp.traj.eq_motion.classical...
                        (obj, pot_mats, frc_mats, D, U, m, psi);
            
            % Add non-classical contributions to adiabatic force
            frc_adi = get_frc_qt_adi_and_mom_dif(obj, pot_mats, frc_mats, D, U, frc_adi, psi);            
        end
        
        % Calculate non-classical contributions to adiabatic force
        function [frc_qt_adi, mom_diff] = get_frc_qt_adi_and_mom_dif(obj, pot_mats, frc_mats, D, U, frc, psi) 
            global space
            frc_qt_adi = cell(space.n_dim,1);
            mom_diff   = cell(space.n_dim,1);
            
            m=1;
            n=2;
            
            nac_mn = ham.nac_mn (pot_mats,frc_mats,U,D,m,n);
            
            coherence = conj ( psi{n} ) .* psi{m};
            
            a = real(coherence);
            b = imag(coherence);
            
            E_diff(:,1) = D(m,:) - D(n,:);
            
            for d = 1:space.n_dim
                frc_qt_adi{d}   = frc{d} + 2 * E_diff .* nac_mn{d} .* a;
                mom_diff{d}     = 2 * nac_mn{d} .* b;
            end  
            
        end
        
        % Calculate difference between "normal" and kinetic momentum
        function mom_diff = get_mom_diff (obj,pot_mats,frc_mats,D,U,psi)
            global space
            mom_diff   = cell(space.n_dim,1);
            
            m=1;
            n=2;
            
            nac_mn = ham.nac_mn (pot_mats,frc_mats,U,D,m,n);
            
            coherence = conj ( psi{n} ) .* psi{m};
            b = imag(coherence);
            
            for d = 1:space.n_dim
                mom_diff{d}     = 2 * nac_mn{d} .* b;
            end  
        end
        
    end
end



    
