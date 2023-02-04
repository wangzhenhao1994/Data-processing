%--------------------------------------------------------------------------
%
% Compute force along m-th potential energy surface 
% from diabatic force matrices and eigenvector matrices
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function frc_adi = frc_adi (pot_mats,frc_mats,U,m)
global space hamilt

n_traj = length(U(1,1,:));

if n_traj>0
    
    frc_adi = cell(space.n_dim,1);
    
    switch hamilt.coupling.n_eqs
        case 1
            frc_adi = frc_mats;
        case 2
            V_11(:,1) = pot_mats(1,1,:);
            V_22(:,1) = pot_mats(2,2,:);
            V_12(:,1) = pot_mats(1,2,:);
            
            dlt = (V_11 - V_22)/2;
            rho = sqrt ( dlt.^2 + V_12.^2 );
            
            % Forces along adiabatic potential curves|surfaces
            % Analytic solutions for a 2-state problem
            for d=1:space.n_dim
                F_11_d(:,1) = frc_mats{d}(1,1,:);
                F_22_d(:,1) = frc_mats{d}(2,2,:);
                F_12_d(:,1) = frc_mats{d}(1,2,:);
                
                flt_d = (F_11_d - F_22_d)/2;
                fta_d = (F_11_d + F_22_d)/2;
                fra_d = ( dlt.* flt_d + V_12.*F_12_d ) ./ rho;
                
                if(m==1)
                    frc_adi{d} = fta_d - fra_d; % lower adiabat
                else
                    frc_adi{d} = fta_d + fra_d; % lower adiabat
                end
            end
            
        otherwise
            
            % Calculate and save adiabatic forces
            % Eq. (7) from doi:10.1007/3-540-45054-8_36
            for d=1:space.n_dim
                % Computation of the adiabatic force: F_adi_d = u_m^T * F_d * u_m
                u_m_F_d = reshape( sum( U(:,m,:) .* frc_mats{d} , 1 ) , [hamilt.coupling.n_eqs,n_traj] );
                u_m     = squeeze(U(:,m,:));
                
                frc_adi{d} = zeros(n_traj,1);
                frc_adi{d}(:) = sum( u_m_F_d .* u_m ,1);
            end
    end
end
        
