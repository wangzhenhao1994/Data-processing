%--------------------------------------------------------------------------
%
% "Single switch surface hopping" class definition
%
% Hopping probabilities are evaluated:
% at crossings of diabatic potentials OR
% at minima of gaps between adiabatic potentials
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

classdef sssh < sht.generic & handle
    
    properties (Access = public)
        lz_variant  % Choice of different LZ variants
        
        D_old       % Adiabatic energies (eigenvalues) of the previous step
        D_oldold    % Adiabatic energies (eigenvalues) of the prev prev step
        
        gap_m1      % Energy gaps at previous time step
        gap_m2      % Energy gaps at pre-previous time step
    end
    
    methods (Access = public)
        
        % Constructor: Set default values for class properties
        function obj = sssh(lzv,n,seed)
            
            % Inherit from "generic" super class
            obj = obj@sht.generic(n,seed);
            
            % Setting LZ = Landau Zener variant to be used
            if isempty(lzv)
                obj.lz_variant = 2;      % energy gaps only, no NAC tensors required 
            else
                obj.lz_variant = lzv;
            end
            obj.consec_hops = true;      % only hops to neighbouring levels are allowed
            
            % Validate input: Only six variants implemented
            if obj.lz_variant>6
                prt.error ('Wrong choice of Landau-Zener variant. Only six variants implemented')
            end
                
            % Setting text strings
            obj.string0 = ['sssh' int2str(obj.lz_variant)];    
            obj.string4 = ['Single switch surface hopping: LZ variant ' int2str(obj.lz_variant)];
        end
        
        function save_previous (obj) 
            
            % Inherit from initialization of superclass
            save_previous@sht.generic ( obj);
            
            obj.D_oldold  = obj.D_old;
            obj.D_old     = obj.D_new;
        end
        
        function probable = prob_hop (obj,m,n,ind_m)
            global space time
            
            % Preallocate
            probable = zeros(size(ind_m));
            
            % Adiabatic energy gap
            if ~isempty(obj.D_oldold)
                gap        = abs( obj.D_new(m,ind_m)   -obj.D_new(n,ind_m)    )';
                gap_old    = abs( obj.D_old(m,ind_m)   -obj.D_old(n,ind_m)    )';
                gap_oldold = abs( obj.D_oldold(m,ind_m)-obj.D_oldold(n,ind_m) )';
                
                % Single switch citerion: Hopping can only occur at critical phase space points
                % Detect sign change of first derivative: from negative to positive
                ind_c = find ( gap_oldold > gap_old & gap > gap_old );
                
                if(~isempty(ind_c))
                    
                    pot_mats = obj.pot_mat(:,:,ind_m(ind_c));
                    
                    mom      = cell(space.n_dim,1);
                    frc_mats = cell(space.n_dim,1);
                    for d = 1:space.n_dim
                        mom{d}      = obj.mom{d}(ind_m(ind_c));
                        frc_mats{d} = obj.frc_mat{d}(:,:,ind_m(ind_c));
                    end
                    
                    % https://en.wikipedia.org/wiki/Finite_difference_coefficient
                    gap_d2  = (gap(ind_c) - 2* gap_old(ind_c) + gap_oldold(ind_c)) / time.steps.s_delta^2;
                    
                    probable (ind_c) = lz_formula(obj, mom,  pot_mats, frc_mats, ...
                                            obj.U_new(:,:,ind_m(ind_c)), obj.D_new(:,ind_m(ind_c)), ...
                                            gap(ind_c), gap_d2, m, n);
                end
            end
        end
      
        % see separate files for the following public methods
        probable = lz_formula(obj, mom, dia_pot_mats, dia_frc_mats, U, D, gap, gap_d2, m, n)
        
    end
        
end
    
