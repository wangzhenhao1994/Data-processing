%--------------------------------------------------------------------------
%
% Initial state of a quantum state "ket" object
% =============================================
%
% Create initial and thermal density matrix for LVNE propagation.
% Note that the matrices are vectorized: default is columnwise.
%
% Note that quantum states start counting from zero!
%
% For Boltzmann populations, temperature is in atomic units of 315,923.5 K
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019 - .... Burkhard Schmidt
%               2011        Boris Schaefer-Bung
%               2012        Jeremy Rodriguez, Ulf Lorenz
%
% see the README file for license details.

function init_obj (rho,H0)
global control time

%% When this function is called by qm_abncd
if nargin>1
    
    prt.disp (' ')
    prt.disp ('-------------------------------------------------------------')
    prt.disp (' ')
    
    if ~isfield(time, 'rho')
        prt.error ('No information on initial rho (density matrix) found')
    end
    
    % Size of state vectors
    dim = size(H0,1);
    
    %% Initial density matrix
    mat = zeros(dim);
    switch lower (time.rho.choice)
        case 'pure' % pure state
            ii = time.rho.pure+1;
            mat(ii,ii) = 1;
            prt.disp(['LvNE initial density: pure state = ', int2str(time.rho.pure)])
        case 'cat' % Schroedinger cat state: coherent superposition
            ii = time.rho.cat(1)+1;
            jj = time.rho.cat(2)+1;
            mat(ii,ii) = 1/2;
            mat(jj,jj) = 1/2;
            mat(ii,jj) = 1/2;
            mat(jj,ii) = 1/2;
            prt.disp(['LvNE initial density: cat state = ', int2str(time.rho.cat)])
        case 'mixed' % Incoherent superposition
            ii = time.rho.mixed(1)+1;
            jj = time.rho.mixed(2)+1;
            mat(ii,ii) = 1/2;
            mat(jj,jj) = 1/2;
            prt.disp(['LvNE initial density: mixed state = ', int2str(time.rho.mixed)])
        case 'thermal' % thermal (Boltzmann) distribution
            if (time.rho.temperature==0)
                mat(1,1) = 1;
            else
                boltz = exp(-H0/time.rho.temperature);
                mat = diag(boltz/sum(boltz));
            end
            prt.disp(['LvNE initial density: thermal with kBT = ', num2str(time.rho.temperature)])
            for n=1:dim
                prt.disp([int2str(n-1) ' : ' num2str(real(boltz(n)))])
            end
        otherwise
            prt.error (['Wrong choice of LvNE initial conditions: ' time.rho.choice])
    end
    prt.disp(' ')
    
    % vectorize initial/thermal density matrices: columnwise ordering
    rho.x_initial = mat(:);
    
    
    %% Thermal density matrix as fixpoint of matrix A ==> Equilibrium
    mat = zeros(dim);
    if (control.lvne.temperature==0)
        mat(1,1) = 1;
    else
        boltz = exp(-H0/control.lvne.temperature);
        mat = diag(boltz/sum(boltz));
    end
    prt.disp(['LVNE equilibrium density: thermal with kBT = ', num2str(control.lvne.temperature)])
    for n=1:dim
        prt.disp([int2str(n-1) ' : ' num2str(mat(n,n))])
    end
    prt.disp(' ')
    
    % vectorize initial/thermal density matrices: columnwise ordering
    rho.x_equilib = mat(:);
    
    %% When this function is called by qm_propa
else
    
    load_0 (rho, 1)
    rho.x = rho.x_initial;
    
end


