%--------------------------------------------------------------------------
%
% Matrix (eigen) representations of important operators using
% output from previous qm_bound calculation
%
% This function calculates the eigenenergies, matrix elements of dipole
% and/or polarizability along x and/or y (if applicable) and optionally 
% also matrix elements of the system/bath coupling and/or additional multi-
% plicative operators. These vectors/matrices are added to object state (of  
% class wave) and are written to file "wave_0.mat" in the current directory. 
%
% An optional second parameter can be supplied to define a cut-off; matrix
% elements are set to zero if their absolute values are below the cut-off, 
% which keeps the output files more readable (and could be used for sparse
% representations later on).
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2014 - 2017 Burkhard Schmidt
%               2011 Ulf Lorenz, Boris Schaefer-Bung, Burkhard Schmidt
%               2012 Jeremy Rodriguez, Burkhard Schmidt, Ulf Lorenz
%               
%
% see the README file for license details.

function qm_matrix ( cutoff )

global control hamilt state time

% Initializes general information and sets up log files.
prt.init (mfilename('fullpath'));

% Calculations of matrix elements for wavefunctions only
if ~isa(state,'wave') 
    prt.error ('Calculations of matrix elements for wavefunctions only')
end

% Provide default values for missing input arguments
if nargin==0
    cutoff=0;
end

prt.disp ('***************************************************************')
prt.disp ('Matrix representations of important operators     ')
prt.disp ('***************************************************************')
prt.disp (' ')
prt.disp ('by numerical integration using DVR (quadrature) schemes     ')
prt.disp ('https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_matrix ');
prt.disp (' ')

%% Load the input
load_0 (state, 2);

%% Preallocate fields (ham, dip, pol, sbc, amo) of structure 'tise' (if applicable)
state.M_ham = zeros(time.steps.m_number, 1);
prt.disp('Energies: \delta_ij E_i = <i|H|j> ')
prt.disp(' ')

% Dipole interactions
if isfield(hamilt,'dip')
    for p = 1:length(hamilt.dip)
        if ~isempty(hamilt.dip{p}{1,1}.dvr)
            state.M_dip{p} = zeros(time.steps.m_number);
            prt.disp(['Dipole moments: \mu_{p,ij} = <i|\mu_p|j> with p=' int2str(p)])
            prt.disp(' ')
        else
            state.M_dip{p} = [];
        end
    end
end

% Polarizability interactions
if isfield(hamilt,'pol')
    for p = 1:size(hamilt.pol,1)
        for q = p:size(hamilt.pol,2)
            if ~isempty(hamilt.pol{p,q}{1,1}.dvr)
                state.M_pol{p,q} = zeros(time.steps.m_number);
                prt.disp(['Polarizabilities: \alpha_{pq,ij} = <i|\alpha_pq|j> with pq='  int2str(p)  int2str(q)])
                prt.disp(' ')
            else
                state.M_pol{p,q} = [];
            end
        end
    end
end

% System-bath couplings
if isfield(hamilt,'sbc')
    if ~isempty(hamilt.sbc{1,1}.dvr)
        state.M_sbc = zeros(time.steps.m_number);
        prt.disp('System-bath coupling: \chi_{ij} = <i|\chi|j> ')
        prt.disp(' ')
    end
end

% Additional multiplicative operators
if isfield(hamilt, 'amo') && strcmpi(control.observe.types,'amo')
    state.M_amo = cell (length(hamilt.amo),1);
    if ~isempty (hamilt.amo{p})
        for p = 1:length(hamilt.amo)
            if ~isempty(hamilt.amo{p}.dvr)
                state.M_amo{p} = zeros(time.steps.m_number);
                prt.disp(['Additional multiplicative operators: O_{p,ij} = <i|O_p|j> with O being: ' hamilt.amo{p}.label ])
                prt.disp(' ')
            end
        end
    end
end

%% Calculate all matrix elements

% Outer loop: bra-states
for bra_index = 1:time.steps.m_number
	load ( state, bra_index );
	bra_state = state.dvr;
    
    % calculate the eigenenergies
    apply_ham( state, [0 0], 0 ); % provides H|psi> in psi.new
    for m = 1:hamilt.coupling.n_eqs
        state.M_ham(bra_index) = wave.braket(bra_state, state.new);
    end
    
	% Inner loop: ket-states 
	for ket_index = 1:time.steps.m_number
		load( state, ket_index );
		ket_state = state.dvr;

        % dipole moment
        if ~isempty (state.M_dip)
            for p = 1:length(hamilt.dip)
                if ~isempty(hamilt.dip{p}{1,1}.dvr)
                    state.M_dip{p}(bra_index, ket_index) = wave.sandwich(bra_state, hamilt.dip{p}, ket_state);
                end
            end
        end

        % polarization
        if ~isempty (state.M_pol)
            for p = 1:size(hamilt.pol,1)
                for q = p:size(hamilt.pol,2)
                    if ~isempty(hamilt.pol{p,q}{1,1}.dvr)  
                        state.M_pol{p,q}(bra_index, ket_index) = wave.sandwich(bra_state, hamilt.pol{p,q}, ket_state);
                    end
                end
            end
        end
        
        % system-bath coupling
        if ~isempty (state.M_sbc)
            state.M_sbc(bra_index, ket_index) = wave.sandwich(bra_state, hamilt.sbc, ket_state);
        end
        
        % additional multiplicative operators
        if ~isempty(state.M_amo)
            for p = 1:length(hamilt.amo)
                if ~isempty (hamilt.amo{p})
                    if ~isempty(hamilt.amo{p}.dvr)
                        state.M_amo{p}(bra_index, ket_index) = wave.sandwich(bra_state, hamilt.amo{p}, ket_state);
                    end
                end
            end
        end
        
    end
end

%% Optionally truncate matrix elements
if cutoff>0
    if ~isempty (state.M_dip)
        for p = 1:length(state.M_dip)
            state.M_dip{p}(abs(state.M_dip{p}) < cutoff) = 0;
        end
    end
    
    if ~isempty (state.M_pol)
        for p = 1:size(state.M_pol,1)
            for q = p:size(state.M_pol,2)
                state.M_pol{p,q}(abs(state.M_pol{p,q}) < cutoff) = 0;
            end
        end
    end
    
    if ~isempty (state.M_sbc)
        state.M_sbc(abs(state.M_sbc) < cutoff) = 0;
    end
    
    if ~isempty (state.M_amo)
        for p = 1:length(hamilt.amo)
            if ~isempty (hamilt.amo{p})
                state.M_amo{p}(abs(state.M_amo{p}) < cutoff) = 0;
            end
        end
    end
    
end

%% Matrix/vector representations and labels of observables (to be used in qm_abncd)
state.M_lab = cell(length(control.observe.choices),1);
state.M_obs = control.observe.types;

switch lower(control.observe.types)
    
    % Additional multiplicative operators
    case 'amo'
        if isempty (state.M_amo)
            prt.error ('No additional multiplicative operators defined')
        end
        
        % Only observables corresponding to *one* operator
        for len=1:length (control.observe.choices)
            if length(control.observe.choices{len})>1
                prt.error('Combinations of more than one observables not *yet* implemented')
            end
        end

        % If not specified otherwise, observables will be labeled like operators
        if ~isfield (control.observe,'labels')
            for len=1:length (control.observe.choices)
                control.observe.labels{len} = hamilt.amo{len}.label;
            end
        end
            
        % Set labels and observables ==> structure "tise"
        state.M_mat = cell(length(control.observe.choices),1);
        for len=1:length (control.observe.choices)
            state.M_lab{len} = control.observe.labels{len};
            prt.disp (['Observable ' int2str(len) ': Additional multiplicative operators: ' state.M_lab{len}])
            state.M_mat{len} = state.M_amo{control.observe.choices{len}};
        end
        
    % Populations as projectors onto eigenstates
    case 'prj'
        state.M_mat = cell(length(control.observe.choices),1);
        for len=1:length(control.observe.choices)
            state.M_lab{len} = control.observe.labels{len};
            prt.disp (['Observable ' int2str(len) ': Populations of eigenstates: ' state.M_lab{len}])
            state.M_mat{len} = zeros (time.steps.m_number);
            for m=1:length(control.observe.choices{len})
                ii=control.observe.choices{len}(m)+1;
                state.M_mat{len}(ii,ii) = 1;
            end
        end
        
    % Populations from overlaps with eigenstates
    case 'ovl'
        state.M_vec = cell(length(control.observe.choices),1);
        for len=1:length(control.observe.choices)
            state.M_lab{len} = control.observe.labels{len};
            prt.disp (['Observable ' int2str(len) ': Overlaps with eigenstates: ' state.M_lab{len}])
            state.M_vec{len} = zeros (time.steps.m_number,1);
            for m=1:length(control.observe.choices{len})
                ii=control.observe.choices{len}(m)+1;
                state.M_vec{len}(ii,1) = 1;
            end
        end
        
    otherwise
        prt.error('Wrong choice of observable types')
end

%% Save results to file
prt.disp (' ')
save_0(state); 

%% Plot matrices

% Open new figure
figure(11);
clf;
thisplot = vis.styles; % Construct object
show_logo (thisplot)

% Energy level diagram (diagonal matrix)
subplot(2,2,1)
plot(0:time.steps.m_number-1,real(state.M_ham),'o')
set(gca, ...
    'LineWidth',  thisplot.l_thick, ...
    'FontName',   thisplot.f_name, ...
    'FontSize',   thisplot.f_large, ...
    'FontWeight', thisplot.f_heavy)
xlabel('n')
ylabel('E(n)')
title('Energy level diagram')

% System-bath coupling
if ~isempty (state.M_sbc)
    subplot(2,2,2)
    surf(state.M_sbc)
    set(gca, ...
        'LineWidth',  thisplot.l_thick, ...
        'FontName',   thisplot.f_name, ...
        'FontSize',   thisplot.f_large, ...
        'FontWeight', thisplot.f_heavy)
    xlabel('n')
    ylabel('m')
    zlabel('\chi(n,m)')
    title('System-bath coupling')
end

% Dipole moments (along x and/or y)
if ~isempty (state.M_dip)
    for p=1:length(state.M_dip)
        if ~isempty (state.M_dip{p})
            subplot(2,2,2+p)
            surf(state.M_dip{p})
            set(gca, ...
                'LineWidth',  thisplot.l_thick, ...
                'FontName',   thisplot.f_name, ...
                'FontSize',   thisplot.f_large, ...
                'FontWeight', thisplot.f_heavy)
            xlabel('n')
            ylabel('m')
            switch p
                case 1
                    zlabel('\mu_x(n,m)')
                    title('Dipole moments (x)')
                case 2
                    zlabel('\mu_y(n,m)')
                    title('Dipole moments (y)')
            end
        end
    end
end


% Output clock/date/time
prt.clock;

end

