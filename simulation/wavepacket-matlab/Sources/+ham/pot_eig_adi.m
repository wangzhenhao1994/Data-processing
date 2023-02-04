%--------------------------------------------------------------------------
% 
% Compute adiabatic potentials, i.e. eigenvalues and 
% eigenvector matrix of diabatic potential matrix
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function [U,D] = pot_eig_adi(pot_mats)
global hamilt time

n_traj = length(pot_mats(1,1,:));
U = zeros ( hamilt.coupling.n_eqs , hamilt.coupling.n_eqs , n_traj );
D = zeros ( hamilt.coupling.n_eqs , n_traj );

switch hamilt.coupling.n_eqs
    case 1
        D(1,:) = pot_mats(1,1,:);
    case 2
        V_11(:,1) = pot_mats(1,1,:);
        V_22(:,1) = pot_mats(2,2,:);
        V_12(:,1) = pot_mats(1,2,:);
        
        % Adiabatic potential energy curves|surfaces
        % Analytic solutions for a 2-state problem
        % See e.g. Eq. (4.8) in doi:10.1063/1.1522712
        dlt = (V_11 - V_22)/2;
        eta = (V_11 + V_22)/2;
        rho = sqrt ( dlt.^2 + V_12.^2 );
        
        D(1,:) = eta - rho; % lower adiabat
        D(2,:) = eta + rho; % upper adiabat
        
        % Compute eigenvectors with explicit formulas (not for SSSH, variant 2)
        ind_n = V_12~=0;
        n21   = (V_11 > V_22)*1.0;
        n22   = (V_11 < V_22)*1.0;
        
        b21          = V_11(ind_n)-V_22(ind_n)+2*rho(ind_n);
        b22          = 2*V_12(ind_n);
        n22(ind_n)   = b22 ./ sqrt(b21.^2 + b22.^2);
        n21(ind_n)   = b21 ./ sqrt(b21.^2 + b22.^2);
        
        % Save previous and new eigenvectors
        U(1,1,:) = - n22;
        U(2,1,:) =   n21;
        U(1,2,:) =   n21;
        U(2,2,:) =   n22;
        
    otherwise
        
        % Diagonalize diabatic potential matrices for each trajectory
        for t = 1:n_traj
            % Entries of vector D are eigenvalues
            % Columns of matrix U are right eigenvectors
            [U(:,:,t),D(:,t)] = eig(pot_mats(:,:,t),'vector');
        end
end

time.counter.diag.total = time.counter.diag.total + n_traj;

end

