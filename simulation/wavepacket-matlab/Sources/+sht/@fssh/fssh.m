%--------------------------------------------------------------------------
%
% FSSH = Fewest switching surface hopping algorithm
% 
% Determine transition probability from d/dt |c(t)|^2
% from TDSEs attached to each of the trajectories.
% According to the proof given by Tully, this algorithm
% really results in the lowest number of switches.
% 
% see: J. C. Tully
%      J. Chem. Phys. 93(2), 1061-1071 (1990)
%      DOI:10.1063/1.459170
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef fssh < sht.mssh & handle
    
    methods (Access = public)

        % Constructor: Setting defaults and a text string
        function obj = fssh (n,seed)
            obj = obj@sht.mssh(n,seed);     % Inherit from superclass
            obj.string0 = 'fssh';
            obj.string4 = 'Fewest switches surface hopping';
        end   
                
        % Get probabilities of hopping from state "m" to state "n"
        % Eqs. (14,15) from doi:10.1063/1.459170
        function probable = prob_hop (obj,m,n,ind_m)
            global hamilt space time
            
            % Preallocate
            probable = zeros(size(ind_m));
            
            % Quantum coherence
            coherence = conj ( obj.psi{n}(ind_m) ) .* obj.psi{m}(ind_m);
            
            % Diabatic picture
            if strcmpi(hamilt.coupling.represent,'dia')
                coupling = obj.ham {n,m}(ind_m);
                probable = +2 * imag (coupling .* coherence);
                
            % Adiabatic picture
            elseif strcmpi(hamilt.coupling.represent,'adi')
                
                % Quantities needed to calculate NAC vector
                mom      = cell(space.n_dim,1);
                frc_mats = cell(space.n_dim,1);
                for d = 1:space.n_dim
                    mom{d}      = obj.mom{d}(ind_m);
                    frc_mats{d} = obj.frc_mat{d}(:,:,ind_m);
                end
                U = obj.U_new(:,:,ind_m);
                D = obj.D_new(:,ind_m);
                pot_mats = obj.pot_mat(:,:,ind_m);
                
                % Calculate NAC vector
                nac_mn = ham.nac_mn (pot_mats,frc_mats,U,D,m,n);
                
                % Calculate coupling (see FSSH formula)
                coupling = zeros (size(probable));
                for d = 1:space.n_dim
                    coupling = coupling + nac_mn{d} .* mom{d} / space.dof{d}.mass;
                end
                
                % Compute FSSH formula
                probable = 2 * real (coupling .* coherence);
                
            end
            probable = time.steps.s_delta * probable ./ abs(obj.psi{m}(ind_m)).^2;
            
        end
        
    end
end

