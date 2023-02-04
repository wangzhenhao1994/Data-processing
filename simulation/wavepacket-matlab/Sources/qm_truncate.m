%------------------------------------------------------------------------------
%
% Truncate/reduce A, B, N, C and initial/equililibrium state vectors
%
% Upon balancing transformation, the obtained Hankel singular values (HSVs)
% indicate controllability and observability, and the resulting modes are 
% ordered accordingly. To reduce dimensionality, there are two methods: 
%
% Method = 't' (truncatation): 
% In simple tuncation one simply eliminates the low HSV modes because
% they are weakly controllable and weakly observable at the same time. 
%
% Method = 's' (singular perturbation):
% The idea of confinement truncation is based on the analogy 
% between large HSV-modes with slow dof's and low HSV-modes
% with fast dof's. Then an averaging principle ( based on singular 
% perturbation theory) can be used to derive equations of motion 
% for the former one, where the latter one is confined to its average, 
% or rather, the t --> oo limit.  
%
% A = A11 - A12 A22 \ A21
% N = A11 - N12 A22 \ A21
% C = C1  - C2  A22 \ A21
%
% From both an execution time and numerical accuracy standpoint, it is better
% to use Matlab's matrix division operator mldivide (A22\...) than inv(A22)
%
% Already balanced A, B, N, C matrices and initial/equilibrium state vectors 
% are read from "name_bal_0.mat" where name typically is 'ket' or 'rho'. 
% In the end, the truncated matrices and vectors are written to either 
% "name_t<n>_0.mat" (method='t') for (simple) truncation or 
% "name_s<n>_0.mat" (method='s') for singular perturbation theory.
% Input parameter 'dim' specifies the dimension <n> of the truncated system.
%
% C. Hartmann, B. Schaefer-Bung, A. Thöns-Zueva
% SIAM J. Control Optim., 51(3), 2356-2378 (2013)
% doi:10.1137/100796844
% 
% C. Hartmann, V. Vulcanov, Ch. Schütte
% Multiscale Model. Simul., 8(4), 1348-1367 (2010)
% doi:10.1137/080732717
% 
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Boris Schaefer-Bung
% Copyright (C) 2012-16 Burkhard Schmidt
%
% see the README file for license details.
% 

function qm_truncate (method, dim)

global state

% Initializes general information and sets up log files.
prt.init (mfilename('fullpath'));

% Truncation for state vectors and density matrices only
if ~isa(state,'ket') 
    prt.error ('Truncation for "ket", "rho" only')
end

% Load full matrices (already balanced) from data file
state.save_suffix = 'bal';
load_0 (state,true);

prt.disp (' ')
prt.disp ('-------------------------------------------------------------')
prt.disp (' Truncating A, B, N, C matrices; also x_i, x_e state vectors ')
prt.disp (' http://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_truncate');
prt.disp ('-------------------------------------------------------------')
prt.disp (' ')
prt.disp(['truncating to ' int2str(dim) ' modes:'])

nbal=size(state.A,1);

switch lower(method)
    case 's' % singular perturbation for A, N, C
        
        prt.disp('singular perturbation approximation')
        A = state.A;
        state.A = A(1:dim,1:dim) ...                          % A11
                - A(1:dim,dim+1:nbal) ...                     % A12
                *(A(dim+1:nbal,dim+1:nbal) ...                % A22^(-1)
                \ A(dim+1:nbal,1:dim));                       % A21
            
        N = state.N;
        for d=1:length(N)
            state.N{d} = N{d}(1:dim,1:dim) ...                % N11
                       - N{d}(1:dim,dim+1:nbal) ...           % N12
                       *(A(dim+1:nbal,dim+1:nbal) ...         % A22^(-1)
                       \ A(dim+1:nbal,1:dim));                % A21
        end
        
        C = state.C;
        for d=1:length(C)
            state.C{d} = C{d}(1:dim) ...                      % C1
                       - C{d}(dim+1:nbal) ...                 % C2
                       *(A(dim+1:nbal,dim+1:nbal) ...         % A22^(-1)
                       \ A(dim+1:nbal,1:dim));                % A21
        end
        
        state.save_suffix = ['s' int2str(dim)];                % file name suffix
        state.title=[state.title 'perturb. reduction to ']; % simulation title
        
    case 't' % truncation  for A, N, C
        
        prt.disp('simple truncation')
        state.A=state.A(1:dim,1:dim);
        for d=1:length(state.N)
            state.N{d}=state.N{d}(1:dim,1:dim);
        end
        for d=1:length(state.C)
            state.C{d}=state.C{d}(1:dim);
            state.Q{d}=state.Q{d};
        end

        state.save_suffix = ['t' int2str(dim)];                % file name suffix
        state.title=[state.title 'simple truncation to '];  % simulation title
        
    otherwise
        util_error (['Invalid choice of truncation method : ' method])
end

check_stable (state, 'truncated system matrix A')

% truncate/reduce B matrices
for d=1:length(state.B)
    state.B{d}=state.B{d}(1:dim,:);
end

% truncate transformation matrices and state vectors (initial and equilibrium) 
state.S=state.S(1:dim,:); % needed for qm_correlate ==> oct.reconstruct
state.T=state.T(:,1:dim); % needed for qm_correlate ==> oct.reconstruct
state.x_initial=state.x_initial(1:dim,:);
state.x_equilib=state.x_equilib(1:dim,:);

% Save truncated/reduced matrices to data file
save_0 (state)

% Plot spectrum of A 
spectrum_A (state, mfilename)

% Output clock/date/time
prt.clock;

end

