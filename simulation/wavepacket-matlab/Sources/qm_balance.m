%------------------------------------------------------------------------------
%
% Balancing transformation of controllability and observability Gramians:
% A, B, N, C matrices and initial/equilibrium state vectors are read from
% data file where the filename typically is 'rho_0.mat'. In the end, the
% the balanced matrices and vectors are written to file 'rho_bal_0.mat'.
%
% For more details on the background of the balancing transformation see
% B. Schaefer-Bung, C. Hartmann, B. Schmidt, and Ch. Schuette
% J. Chem. Phys. 135, 014112 (2011)
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Boris Schaefer-Bung
%               2013-16 Burkhard Schmidt

function qm_balance

global reduce state
 
% Initializes general information and sets up log files.
prt.init (mfilename('fullpath'));

prt.disp (' ')
prt.disp ('------------------------------------------------------------------------------')
prt.disp (' Balancing transformation of A, B, N, C matrices; also x_i, x_e state vectors ')
prt.disp (' http://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_balance')
prt.disp ('------------------------------------------------------------------------------')
prt.disp (' ')

% Set all parameters to default values if not specified by the user 
if ~isfield (reduce,'balance')
    reduce.balance=[];
end

if ~isfield (reduce,'glyap')
    reduce.glyap=[];
end

if ~isfield (reduce.balance,'method')
    reduce.balance.method='iter';
end
prt.disp(['Method for solving generalized Lyapunov equation: ' reduce.balance.method])

if ~isfield (reduce.balance,'transform')
    reduce.balance.transform = 'srbt';
end
prt.disp(['Method of balancing transformation: ' reduce.balance.transform])

if ~isfield(reduce.balance,'BN_scale')
    reduce.balance.BN_scale = 1;
end
prt.disp(['Scaling B-vector, N-matrix: ' num2str(reduce.balance.BN_scale)])
prt.disp(' ')

% Balancing transformation for state vectors and density matrices only
if ~isa(state,'ket') 
    prt.error ('Balancing transformation for "ket", "rho" only')
end

% Load full matrices from data file  
load_0 (state,true)
n = size(state.A,1);
prt.disp (['Dimensionality of original system: ' int2str(n)])

% Retrieve ABNCD information from object "state" and downscale
% B and N such that generalized Lyapunov equations are solvable
A = state.A;

B = cell(length(state.B),1);
for d=1:length(state.B)
    B{d} = state.B{d}/reduce.balance.BN_scale;
end

N = cell(length(state.N),1);
for d=1:length(state.N)
    N{d} = state.N{d}/reduce.balance.BN_scale;
end

C = cell(length(state.C),1);
for d=1:length(state.C)
    C{d}=state.C{d};
end
prt.disp ('   ')

% Check stability: real part of all eigenvalues should be negative
check_stable (state, 'original system matrix A')
prt.disp (' ')

%% Two choices of stabilizing A matrix
prt.disp ('Stabilizing A matrix')
if ~isfield(reduce.balance,'A_stable')
    reduce.balance.A_stable = 'ssu';
end
switch lower(reduce.balance.A_stable)
    case 'ssu' 
        % splitting stable from unstable part i.e. where Re(eig(A)) >= 0
        % consider only the orthogonal complement of zero eigenspace
        % generalization: split all eigenvalues with real part below
        % specified threshold
        
        prt.disp ('Method ''ssu'' : splitting stable from unstable part')     
        if ~isfield (reduce.balance,'A_split')
            reduce.balance.A_split = 1;
        end
        n1 = reduce.balance.A_split; 
        prt.disp (['Dimension of unstable part : ',int2str(n1)])
        prt.disp (' ')
        
        % Calculate all eigenvalues and eigenvectors (columns of U) of A.
        % Note that in general (T>0, LvNE), A is complex and non-symmetric
        [U,D] = eig(full(A));
        [~, ind] = sort(real(diag(D)), 'descend');
        U = U(:,ind);
        
        % Transformation and splitting of A, B, N, C into eigenvector basis
        A  = U\A*U;             % transformation; should be diagonal matrix
        A  = math.sparse(A);    % reinforce sparsity of A
        Au = A(n1+1:n,   1:n1); % unstable part
        A  = A(n1+1:n,n1+1:n ); % stable part
        
        for d=1:length(B)
            B{d} = U\B{d};       % transformation
            B{d} = B{d}(n1+1:n); % stable part
        end
         
        Nu = cell(length(N),1);
        for d=1:length(N) 
            N{d}  = U\N{d}*U;             % transformation
            N{d}  = math.sparse(N{d});    % reinforce sparsity of A
            Nu{d} = N{d}(n1+1:n,   1:n1); % unstable part
            N{d}  = N{d}(n1+1:n,n1+1:n ); % stable part
        end
        
        for d=1:length(C)
            C{d} = C{d}*U;       % transformation
            C{d} = C{d}(n1+1:n); % stable part 
        end
        
        % Rho coupling terms act as additional control field
        % Does it really make things any better?
        if ~isfield (reduce.balance,'acf_couple')
           reduce.balance.acf_couple=0; 
        end
        if reduce.balance.acf_couple
            nBg = length(B);
            B{nBg+1} = Au;
            for d=1:length(N) 
                B{nBg+1+d} = Nu{d};
            end
        end
        
    case 'evs' 
        prt.disp('Method ''evs'': eigenvalue shift of A')
        
        % value of shift
        if ~isfield(reduce.balance,'A_shift')
            reduce.balance.A_shift = 1e-6;
        end
        prt.disp (['Value of shift : ' num2str(reduce.balance.A_shift)])
        prt.disp (' ')
        
        % shift diagonal values of A-matrix
        A = A - reduce.balance.A_shift * speye(n);
                
    otherwise
        prt.error (['Wrong choice of stabilization method : ' reduce.balance.A_stable])
        
end

% Check stability of the splitted (SSU) or shifted (EVS) system
switch lower(reduce.balance.A_stable)
    case 'ssu' 
        check_stable (state, 'splitted system matrix A')
        % check_stable ('splitted system matrix A+N*N^T', A, N)
    case 'evs' 
        check_stable (state, 'shifted system matrix A')
        % check_stable ('shifted system matrix A+N*N^T', A, N)
end
prt.disp (' ')
  
%% Controllability Gramian from generalized Lyapunov equation
prt.disp('------------------------------------');
prt.disp('Solve generalized Lyapunov equation ');
prt.disp('to get controllability Gramian W_C  ');
prt.disp('A*WC + WC*Ap + N*WC*Np + B*Bp = 0   ');
switch lower (reduce.balance.method)
    case 'iter'
        WC = math.glyap1(A,B,N,false,reduce.glyap);
    case 'bicg'
        WC = math.glyap2(A,B,N,false,reduce.glyap);
    otherwise
        prt.error (['Wrong choice of GLYAP solver method : ' reduce.balance.method])
end

%% Observability Gramian from generalized Lyapunov equation
prt.disp('------------------------------------');
prt.disp('Solve generalized Lyapunov equation ');
prt.disp('to get observability Gramian W_O    ');
prt.disp('Ap*WO + WO*A + Np*WO*N + Cp*C = 0   ');
switch lower (reduce.balance.method)
    case 'iter'
        WO = math.glyap1(A,C,N,true,reduce.glyap);
    case 'bicg'
        WO = math.glyap2(A,C,N,true,reduce.glyap);
end

%% Determine method of balancing transform
switch lower (reduce.balance.transform)
    case 'srbt' % Square Root Balanced Truncation  
        [S, T, Sigma] = oct.srbt(WC, WO); 
    case 'mrmr' % Minima Realization and Model Reduction
        [S, T, Sigma] = oct.mrmr(WC, WO);
    otherwise
        prt.error (['Wrong choice of balancing transform method : ' reduce.balance.transform])
end

% Plot details of balancing transformation
vis.balance (Sigma,WC,WO,S,T)

% Add extra blocks to S, T matrices to compensate for the splitting
% Apparently, this leaves T to be right-inverse of S, but not left-inverse 
if strcmpi(reduce.balance.A_stable,'ssu')
    dimt=size(T);
    S=[eye(n1) zeros(n1,dimt(1)); zeros(dimt(2),n1) S];
    S=S/U;
    T=[eye(n1) zeros(n1,dimt(2)); zeros(dimt(1),n1) T];
    T=U*T;
end

%% Carry out balancing transformation
% Check whether T is the (pseudo-)inverse of S 
math.pseudo(S,T)

% Transform A,B,N,C
state.A=S*state.A*T;
for d=1:length(state.B)
    state.B{d}=S*state.B{d};
end
for d=1:length(state.N)
    state.N{d}=S*state.N{d}*T;
end
for d=1:length(state.C)
    state.C{d}=state.C{d}*T;
    state.Q{d}=state.Q{d};
end

% State vectors (x, transformed) 
state.x_initial=S*state.x_initial;
state.x_equilib=S*state.x_equilib;

% Add S and T matrices to object state
state.S=S;  
state.T=T;

% Save balanced matrices etc to data file
state.save_suffix = 'bal';
save_0(state);

% Plot spectrum of A after(!) balancing
spectrum_A (state, mfilename)

% Output clock/date/time
prt.clock;

end
