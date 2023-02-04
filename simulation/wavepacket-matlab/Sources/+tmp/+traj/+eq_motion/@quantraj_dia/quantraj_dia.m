%--------------------------------------------------------------------------
%
% Quantum trajectories in diabatic representation
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

classdef quantraj_dia < tmp.traj.eq_motion.classical & handle
    
    properties (Access = public)
    end    

    methods (Access = public)
        
        % Construct object
        function obj = quantraj_dia
            prt.error('Not yet completely debugged - Not recommended for use');
        end
        
        function solve_tdse(obj, state, step_size, D, U, U_old)
            tdse_adi(state, step_size, D, U, U_old)
        end
        
        % Propagate
        function frc_dia = eq_frc (obj, pot_mats, frc_mats, D, U, m, psi) 
            
            frc_adi = eq_frc@tmp.traj.eq_motion.classical...
                        (obj, pot_mats, frc_mats, D, U, m, psi);
            
            frc_dia = get_frc_qt_dia(obj, frc_mats, U, frc_adi, psi) ;            
        end
        
        function [frc_qt_dia] = get_frc_qt_dia(obj, frc_mats, U, frc, psi) 
            global space
            frc_qt_dia = cell(space.n_dim,1);
            
            m=1;
            n=2;
            
            psi_dia = get_psi_dia (obj,psi,U);
            
            coherence = conj ( psi_dia{n} ) .* psi_dia{m};
            
            a = real(coherence);
            
            for d = 1:space.n_dim
                frc_coupling_dia(:,1) = frc_mats{d}(m,n,:);
                frc_qt_dia{d}   = frc{d} + 2 .* a .* frc_coupling_dia;
            end  
            
        end
        
        function psi_dia = get_psi_dia (obj,psi_adi,U)
            global hamilt
            
            psi_dia = cell(hamilt.coupling.n_eqs,1);
            
            % Computation of the exponential matrix
            for m=1:hamilt.coupling.n_eqs
                psi_dia{m} = zeros ( size(psi_adi{m}) );
                for n=1:hamilt.coupling.n_eqs
                    U_mn(:,1) = U(m,n,:);
                    psi_dia{m} = psi_dia{m} + U_mn .* psi_adi{n};
                end
            end
        end
        
    end
end



    
