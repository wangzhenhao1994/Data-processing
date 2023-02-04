%--------------------------------------------------------------------------
%
% Compute first order non-adiabatic coupling vector between 
% adiatic states m and n from diabatic potential matrices, 
% eigenvectors and eigenvalues thereof, and force matrices.
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function nac_mn = nac_mn (pot_mats,frc_mats,U,D,m,n)
global space hamilt

n_traj = length(U(1,1,:));

if n_traj>0
    
    nac_mn = cell(space.n_dim,1);
    
    switch hamilt.coupling.n_eqs
        case 1
            for d=1:space.n_dim
                nac_mn{d} = zeros(n_traj,1);
            end
        case 2
            V_11 = squeeze(pot_mats(1,1,:));
            V_22 = squeeze(pot_mats(2,2,:));
            V_12 = squeeze(pot_mats(1,2,:));
            
            F_11 = cell(space.n_dim,1);
            F_22 = cell(space.n_dim,1);
            F_12 = cell(space.n_dim,1);
            
            for d=1:space.n_dim
                F_11{d} = squeeze(frc_mats{d}(1,1,:));
                F_22{d} = squeeze(frc_mats{d}(2,2,:));
                F_12{d} = squeeze(frc_mats{d}(1,2,:));
            end
            
            % First order non-adiabatic coupling vectors
            % Analytic solutions for a 2-state problem
            % See e.g. Eq. (4.9) in doi:10.1063/1.1522712
            dlt = (V_11 - V_22)/2;
            rho2 = dlt.^2 + V_12.^2;
            for d=1:space.n_dim
                nac_mn{d} = zeros(n_traj,1);
                
                flt = (F_11{d} - F_22{d})/2;
                C_12_d = ( - dlt.*F_12{d}/2 + V_12.*flt/2 ) ./ rho2;
                
                nac_mn{d}(:) = (-1)^(m>n) * C_12_d;
            end
            
        otherwise
            
            % Calculate non-adiabatic couplings, using also the diabatic(!) forces
            for d=1:space.n_dim
                
                nac_mn{d} = zeros(n_traj,1);
                
                % Computation of u_m^T * F_d * u_n
                F_d = frc_mats{d}(:,:,:);
                u_m_F_d = reshape( sum( U(:,m,:) .* F_d , 1 ) , [hamilt.coupling.n_eqs,n_traj] );
                u_n = squeeze(U(:,n,:));
                u_m_F_d_u_n = sum( u_m_F_d .* u_n ,1);
                
                % Computation of the NAC-vector
                C_mn_d = u_m_F_d_u_n ./ (D(m,:) - D(n,:));
                nac_mn{d}(:) = + C_mn_d ;
                
            end
    end
end
