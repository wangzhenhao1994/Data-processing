%--------------------------------------------------------------------------
%
% Initial state of a quantum state "ket" object
% =============================================
%
% Create initial state vector for TDSE propagation 
%
% Note that quantum states start counting from zero!
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

function init_obj (ket,H0)
global time

prt.disp (' ')
prt.disp ('-------------------------------------------------------------')
prt.disp (' ')

%% When this function is called by qm_abncd
if nargin>1
    
    if ~isfield(time, 'ket')
        prt.error ('No information on initial ket (state vector) found')
    end
    
    % Size of state vectors
    dim = size(H0,1);
    
    % Initial state vector
    vec = zeros(dim,1);
    switch lower (time.ket.choice)
        case 'pure' % pure state
            ii = time.ket.pure+1;
            vec(ii) = 1;
            prt.disp(['TDSE initial state: pure state = ', int2str(time.ket.pure)])
        case 'cat' % Schroedinger cat state
            ii = time.ket.cat(1)+1;
            jj = time.ket.cat(2)+1;
            vec(ii) = 1/sqrt(2);
            vec(jj) = 1/sqrt(2);
            prt.disp(['TDSE initial state: cat state = ', int2str(time.ket.cat)])
        otherwise
            prt.error (['Wrong choice of TDSE initial conditions: ' time.ket.choice])
    end
    ket.x_initial = vec;
    
    % Equilibrium state vector
    ket.x_equilib = zeros(dim,1);
    ket.x_equilib(1) = 1;
    prt.disp(['TDSE "equilibrium" : ground state'])
    prt.disp(' ')
    
    %% When this function is called by qm_propa
else
    
    load_0 (ket, 1)
    ket.x = ket.x_initial;
    
end


