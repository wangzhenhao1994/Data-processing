%--------------------------------------------------------------------------
%
% Provide input for bilinear control problems
%
% Matrix   *A*  is created from energies (and coupling to bath for LvNE)
% Vectors  *B* are created from transition dipole moments
% Matrices *N* are created from transition dipole moments
% Vectors  *C* are created from (linear) observables (LvNE)
% Matrices *D* are created from (quadratic) observables (TDSE)
% Vector *x_initial* is the initial state
% Vector *y_initial* is the corresponding output
% Vector *x_equilib* is the equilibrium state (fix point of A)
% Vector *y_equilib* is the corresponding output
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2014-17 Burkhard Schmidt
%               2012 Boris Schaefer-Bung, Burkhard Schmidt, 
%                    Ulf Lorenz, Jeremy Rodriguez
%
% see the README file for license details.
%

function qm_abncd

global control state

% Initializes general information and sets up log files.
prt.init (mfilename('fullpath'));

% Load the energies, dipole, and (system-bath) coupling matrices.
% Typically those come from previous run of qm_bound and qm_matrix.
% Default is to load from file 'wave_0.mat' in working directory
inp = wave;
if ~isempty(state.save_file)
    inp.save_file = state.save_file;
end
if ~isempty(state.save_dir)
    inp.save_dir = state.save_dir;
end
load_0 (inp,true);
dim=size(inp.M_ham, 1);

switch lower(class(state))
    case 'rho'  
        % Map matrices onto vectors using the so-called tetradic (Liouville) 
        % convention introduced by Mukamel's group
        %
        % We use columnwise ordering of the matrix elements, e.g. for
        % the density rho = (rho_00, rho_10, ..., rho_01, rho_11, ...)
        % Note that this is the standard ordering of Fortran/Matlab anyway
        
        % Initial and equilibrium state vectors/densities for propagation
        init_obj (state, inp.M_ham)

        % set temperature; generate title string for plots
        if ~isfield(control.lvne,'temperature')
            control.lvne.temperature = 0;
        end
        prt.disp (['temperature * k_B = ' num2str(control.lvne.temperature)]) 
        prt.disp (' ') 
        kBT = control.lvne.temperature;
        state.title = ['LvNE: \Theta =' num2str(control.lvne.temperature) ', '];
  
        % Construct omega-matrix (Bohr frequencies)
        w=zeros(dim);
        for k=1:dim
            for l=1:dim
                w(k,l)=inp.M_ham(k)-inp.M_ham(l);
                w = real(w); % to be done: delete?
            end
        end
        
        % Set up matrix A; coherent evolution from omega-matrix
        state.A = diag( reshape(-1i*w, dim^2, 1) );

        % Lindblad operators for population relaxation
        if ~isfield(control,'relax')
            control.relax = [];
        end
        if ~isfield(control.relax,'model')
            control.relax.model = 'none';
        end
        
        % Construct Gamma matrix
        % convention used throughout: Gam(k,l) = Gamma_{k <- l}
        Gam = zeros(dim);
        
        switch lower(control.relax.model)
            
            % Relaxation rates from Fermi's golden rule, see Andrianov & Saalfrank 2006
            case 'fermi'
                prt.disp ('Relaxation rates from Fermi''s golden rule')
                prt.disp (['relaxation rate = ' num2str(control.relax.rate)])
                prt.disp (['==> upper state = ' int2str(control.relax.upper)])
                prt.disp (['==> lower state = ' int2str(control.relax.lower)])
                state.title = [state.title '\Gamma =' num2str(control.relax.rate)];
                state.title = [state.title ' (' int2str(control.relax.upper)];
                state.title = [state.title '->' int2str(control.relax.lower) '), '] ;

                ratios = abs( triu(inp.M_sbc,1) ).^2 / abs( inp.M_sbc(control.relax.lower+1,control.relax.upper+1) )^2;
                
                if kBT == 0 % for zero temperature: only downward relaxation
                    for k=1:dim % upper right triangle of Gamma matrix
                        for l=k+1:dim % l \geq k
                            Gam(k,l)=ratios(k,l) ...
                                *w(control.relax.upper+1,control.relax.lower+1) / w(l,k) ...
                                * control.relax.rate;
                        end
                    end
                    
                else % for finite temperature
                    for k=1:dim % upper right triangle: downward relaxation
                        for l=k+1:dim % l \geq k
                            Gam(k,l) = ratios(k,l) ...
                                * w(control.relax.upper+1,control.relax.lower+1) / w(l,k)...
                                * (1-exp(-w(control.relax.upper+1,control.relax.lower+1)/kBT)) ...
                                / (1-exp(-w(l,k)/kBT)) ...
                                * control.relax.rate;
                            % upward transitions from detailed balance, see Eq. (4)
                            Gam(l,k) = Gam(k,l) * exp(-w(l,k)/kBT);
                        end
                    end
                end
                
            % Relaxation rates from Einstein's spontaneous emission
            case 'einstein'
                prt.disp ('Relaxation rates from Einstein''s spontaneous emission')
                prt.disp (['relaxation rate = ' num2str(control.relax.rate)])
                prt.disp (['==> upper state = ' int2str(control.relax.upper)])
                prt.disp (['==> lower state = ' int2str(control.relax.lower)])
                state.title = [state.title '\Gamma =' num2str(control.relax.rate)];
                state.title = [state.title ' (' int2str(control.relax.upper)];
                state.title = [state.title '->' int2str(control.relax.lower) '), '] ;
                
                ratios = abs( triu(inp.M_dip{1},1) ).^2 / abs( inp.M_dip{1}(control.relax.lower+1,control.relax.upper+1) )^2;
                
                for k=1:dim % upper right triangle: downward relaxation
                    for l=k+1:dim % l \geq k
                        Gam(k,l) = ratios(k,l) ...
                            * (w(l,k) / w(control.relax.upper+1,control.relax.lower+1))^3 ...
                            * control.relax.rate;
                        % upward transitions from detailed balance, see Eq. (4)
                        Gam(l,k) = Gam(k,l) * exp(-w(l,k)/kBT);
                    end
                end
                
%             % Constant relaxation rates (for all downward transitions)
            case 'constant'
                prt.disp ('Constant relaxation rates')
                prt.disp (['relaxation rate = ' num2str(control.relax.rate)])
                state.title = [state.title '\Gamma =' num2str(control.relax.rate)];
                state.title = [state.title ' (const.), '] ;

                for k=1:dim % upper right triangle: downward relaxation
                    for l=k+1:dim % l \geq k
                        Gam(k,l)=control.relax.rate;
                        if kBT>0 % upward transitions from detailed balance
                            Gam(l,k) = Gam(k,l) * exp(-w(l,k)/kBT);
                        end
                    end
                end
                
            % Relaxation rates from data file
            case 'datafile'
                prt.disp ('Reading relaxation rates from file')
                state.title = [state.title '\Gamma from file, '] ;

                data = load ('relax');
                Gam = data.relax;
                
            case 'none'
                prt.disp ('Not accounting for relaxation rates')

            otherwise
                prt.error ('Invalid choice of relaxation model')
                
        end
        
        if ~strcmpi(control.relax.model,'none')
            
            % Construct total dephasing rate, i. e. gamma matrix from Eq. (5)
            GamSummed = sum(Gam, 1);
            gamma.r = zeros(dim);
            for k=1:dim
                for l=1:dim
                    gamma.r(k,l)=0.5*(GamSummed(k)+GamSummed(l));
                end
            end
            
            % population gain, similar to Eq. (A2), but for columnwise order
            for l=0:dim-1
                ll = 1 + l*dim + l;
                for k=0:dim-1
                    kk = 1 + k*dim + k;
                    state.A(ll, kk) = state.A(ll, kk) + Gam(l+1,k+1);
                end
            end
            
            % Include population loss and total dephasing in matrix A
            state.A = state.A + diag( reshape(-gamma.r, dim^2, 1) );
            
            % Find extrema of gamma.r
            min_gr = abs(gamma.r(1,2)); min_lo=1-1; min_up=2-1;
            max_gr = abs(gamma.r(1,2)); max_lo=1-1; max_up=2-1;
            for k=1:dim % upper right triangle with diag: downward relaxation
                for l=k:dim % l \gt k
                    abs_gr = abs(gamma.r(k,l));
                    if abs_gr < min_gr
                        min_gr = abs_gr;
                        min_lo = k-1;
                        min_up = l-1;
                    end
                    if abs_gr > max_gr
                        max_gr = abs_gr;
                        max_lo = k-1;
                        max_up = l-1;
                    end
                end
            end

            prt.disp (['min. relax. dephasing = ' num2str(min_gr)])
            prt.disp (['==> upper state       = ' int2str(min_up)])
            prt.disp (['==> lower state       = ' int2str(min_lo)])
            prt.disp (['==> Bohr frequency    = ' num2str(w(min_up+1,min_lo+1))])
            prt.disp (['max. relax. dephasing = ' num2str(max_gr)])
            prt.disp (['==> upper state       = ' int2str(max_up)])
            prt.disp (['==> lower state       = ' int2str(max_lo)])
            prt.disp (['==> Bohr frequency    = ' num2str(w(max_up+1,max_lo+1))])
                 
        end
        prt.disp (' ')      
        
        % Lindblad operator for pure dephasing
        if ~isfield(control,'depha')
            control.depha = [];
        end
        if ~isfield(control.depha,'model')
            control.depha.model = 'none';
        end
        
        switch lower (control.depha.model)
            
            % Dephasing rates from stochastic Gaussian model (quadratic energy gap dependence)
            case 'gauss'
                prt.disp ('Dephasing rates from stochastic Gaussian model')
                prt.disp (['pure dephasing rate = ' num2str(control.depha.rate)])
                prt.disp (['==> upper state     = ' int2str(control.depha.upper)])
                prt.disp (['==> lower state     = ' int2str(control.depha.lower)])
                prt.disp (['==> Bohr frequency  = ' num2str(w(control.depha.upper+1,control.depha.lower+1))])
                state.title = [state.title '\Gamma^* =' num2str(control.depha.rate)];
                state.title = [state.title ' (' int2str(control.depha.upper)];
                state.title = [state.title '->' int2str(control.depha.lower) '), '] ;

                gamma.d = (w/w(control.depha.upper+1,control.depha.lower+1)).^2 * control.depha.rate;
                
                % Constant dephasing rates
            case 'constant'
                prt.disp ('Constant pure dephasing rates')
                prt.disp (['pure dephasing rate = ' num2str(control.depha.rate)])
                state.title = [state.title '\Gamma^* =' num2str(control.depha.rate)];
                state.title = [state.title ' (const.), '] ;

                gamma.d = ones(size(w)) * control.depha.rate;
                
                % Read dephasing rates from data file
            case 'datafile'
                prt.disp ('Reading pure dephasing data from file')
                state.title = [state.title '\Gamma^* from file, '] ;

                data = load ('depha');
                gamma.d = data.depha;
                
            case 'none'
                prt.disp ('Not accounting for pure dephasing rates')
                
            otherwise
                prt.error ('Invalid choice of pure dephasing model')
                
        end
        
        if ~strcmpi(control.depha.model,'none')
            
            % Contribution to A-matrix
            state.A = state.A + diag( reshape(-gamma.d, dim^2, 1) );
            
            % Find extrema
            min_gd = abs(gamma.d(1,2)); min_lo=1-1; min_up=2-1;
            max_gd = abs(gamma.d(1,2)); max_lo=1-1; max_up=2-1;
            for k=1:dim-1 % upper right triangle: downward relaxation
                for l=k+1:dim % l \geq k
                    abs_gd = abs(gamma.d(k,l));
                    if abs_gd < min_gd
                        min_gd = abs_gd;
                        min_lo = k-1;
                        min_up = l-1;
                    end
                    if abs_gd > max_gd
                        max_gd = abs_gd;
                        max_lo = k-1;
                        max_up = l-1;
                    end
                end
            end
            
            prt.disp (['min. pure dephasing = ' num2str(min_gd)])
            prt.disp (['==> upper state     = ' int2str(min_up)])
            prt.disp (['==> lower state     = ' int2str(min_lo)])
            prt.disp (['==> Bohr frequency  = ' num2str(w(min_up+1,min_lo+1))])
            prt.disp (['max. pure dephasing = ' num2str(max_gd)])
            prt.disp (['==> upper state     = ' int2str(max_up)])
            prt.disp (['==> lower state     = ' int2str(max_lo)])
            prt.disp (['==> Bohr frequency  = ' num2str(w(max_up+1,max_lo+1))])
            
        end
        prt.disp (' ')
        
        % Coupling(s) to control field(s)
        if ~isempty (inp.M_dip)
            for p = 1:length(inp.M_dip)
                if ~isempty (inp.M_dip{p})
                    state.N{p}=zeros(dim^2);
                    
                    % Set up N matrices, similar to  Eq. (A3), but for columnwise order
                    for l=0:dim-1
                        for m=0:dim-1
                            index_lm = 1 + l + m*dim;
                            
                            for k=0:dim-1
                                index_km = 1 + k + m*dim;
                                index_lk = 1 + l + k*dim;
                                
                                state.N{p}(index_lm,index_km) = state.N{p}(index_lm,index_km) + inp.M_dip{p}(l+1,k+1);
                                state.N{p}(index_lm,index_lk) = state.N{p}(index_lm,index_lk) - inp.M_dip{p}(k+1,m+1);
                                
                            end
                            
                        end
                    end
                    
                    % Set up B vectors, similar to Eq. (A5), but for columnwise order
                    state.B{p} = state.N{p}*state.x_equilib;
                    
                end
            end
        end
        
        
        
        % Choice of observables is given in control.observe.targets
        % Set up C vectors: columnwise ordering of matrices
        for len=1:length(control.observe.targets)
            state.label{len} = inp.M_lab{control.observe.targets(len)};
            switch inp.M_obs
                case 'amo'
                    prt.disp (['Observable ' int2str(len) ': Additional multiplicative operators: ' state.label{len}])
                case 'prj'
                    prt.disp (['Observable ' int2str(len) ': Populations as projectors onto eigenstates: ' state.label{len}])
                otherwise
                    prt.error ('Invalid choice of observable for LvNE')
            end
            % Transpose, vectorize, transpose such that C*x gives Tr(O*rho)
            op_mat = inp.M_mat{control.observe.targets(len)};
            op_mat = op_mat.';
            op_vec = op_mat(:);
            state.C {len} = op_vec';
            state.Q {len} = false; % use Re<c|x> as observable

        end
        state.C = state.C'; % make C a row cell vector
        prt.disp (' ')
        
        % if desired: transform full matrix problem
        % columnwise -> diagonal first (cw2df)
        if isfield (control.lvne,'order') && strcmpi(control.lvne.order,'df')
            U=math.cw2df(dim);
            
            state.A = U*state.A*U';
            
            for len=1:length(state.C)
                state.C{len} = state.C{len}*U';
            end
            
            for len=1:length(state.N)
                state.N{len} = U*state.N{len}*U';
                state.B{len} = U*state.B{len};
            end
            
            state.x_initial = U*state.x_initial;
            state.x_equilib = U*state.x_equilib;
        end
        
        % Values of observables for initial and equilibrium state
        for len=1:length(control.observe.targets)
            switch inp.M_obs
                case {'amo','prj'}
                    state.y_equilib(len) = real(state.C{len}*state.x_equilib);
                    state.y_initial(len) = real(state.C{len}*state.x_initial);
            end

        end
        
        % Shift initial state vector w.r.t. its equilibrium state
        state.x_initial = state.x_initial - state.x_equilib;
        
    case 'ket'
        
        % Initial and equilibrium state vectors/densities for propagation
        init_obj (state,inp.M_ham);
        
        % making title string for plots
        state.title = 'TDSE: ' ;
        
        % Set up A; shift the eigenenergies by the ground state energy.
        state.A = -1i*(diag(inp.M_ham)-inp.M_ham(1)*eye(dim));
        
        % Coupling(s) to control field(s)
        if ~isempty (inp.M_dip)
            for p = 1:length(inp.M_dip)
                if ~isempty (inp.M_dip{p})
                    state.N{p}=inp.M_dip{p};
                    state.B{p}=state.N{p}*state.x_equilib;
                end
            end
        end
        
        % Set up C vectors or D matrices
        % Choice of observables is given in control.observe.targets
        for len=1:length(control.observe.targets)
            state.label{len} = inp.M_lab{control.observe.targets(len)};
            switch inp.M_obs
                case 'amo'
                    prt.disp (['Observable ' int2str(len) ': Additional multiplicative operators: ' state.label{len}])
                    state.D{len} = inp.M_mat{control.observe.targets(len)};
                case 'prj'
                    prt.disp (['Observable ' int2str(len) ': Populations as projectors onto eigenstates: ' state.label{len}])
                    state.D{len} = inp.M_mat{control.observe.targets(len)};
                case 'ovl'
                    prt.disp (['Observable ' int2str(len) ': Populations from overlaps with eigenstates: ' state.label{len}])
                    state.C{len} = inp.M_vec{control.observe.targets(len)}';
                    state.Q{len} = true; % use |<c|x>|^2 as observable
                otherwise
                    prt.error ('Invalid choice of observable for TDSE')
            end
        end
        prt.disp (' ')
                
        % Values of observables for initial and equilibrium state
        for len=1:length(control.observe.targets)
            switch inp.M_obs
                case {'amo', 'prj'}
                    state.y_equilib(len) = dot ( state.x_equilib, state.D{len} * state.x_equilib );
                    state.y_initial(len) = dot ( state.x_initial, state.D{len} * state.x_initial );
                case 'ovl'
                    state.y_equilib(len) = abs ( state.C{len}*state.x_equilib )^2;
                    state.y_initial(len) = abs ( state.C{len}*state.x_initial )^2;
            end
        end
        
        % Shift initial state vector w.r.t. its equilibrium state
        state.x_initial = state.x_initial - state.x_equilib;
                
    otherwise
        prt.error (['Wrong choice of class of object "state" : ' class(state)])
        
end

% Labels of the observables
for len=1:length(control.observe.targets)
    state.y_label{len} = inp.M_lab{len};
end


% Sparsity of (complex) matrix A, if density below threshold
density = nnz(state.A)/numel(state.A);
prt.disp (['Density of matrix A: ' num2str(density)])
if density<0.1
    state.A = sparse(state.A);
    prt.disp ('  ==> using sparse storage scheme')
else
    prt.disp ('  ==> using full storage scheme')
end

% Sparsity of (real) matrices N, if density below threshold
for d=1:length(state.N)
    density = nnz(state.N{d})/numel(state.N{d});
    prt.disp (['Density of matrix N_'  int2str(d) ': ' num2str(density)])
    if density<0.1
        state.N{d} = sparse(state.N{d});
        prt.disp ('  ==> using sparse storage scheme')
    else
        prt.disp ('  ==> using full storage scheme')
    end
end
prt.disp(' ')

% Plot spectrum of A 
spectrum_A (state, mfilename)

% Save ABNCD etc matrices to data file where the default is to
% save to file 'ket_0.mat' or 'rho_0.mat' in working directory
save_0 (state)  

% Output clock/date/time
prt.clock;

end

