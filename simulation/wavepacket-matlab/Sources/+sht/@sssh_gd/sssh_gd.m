%--------------------------------------------------------------------------
%
% "Single switch surface hopping" class definition
% using the GD="gradient descent" method to exactly 
% localize (global or local) minima of energy gaps.
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef sssh_gd < sht.sssh & handle
    
    properties (Access = public)
        
        hop_pos     % Position at hopping time 
        hop_mom     % Momentum at hopping time 
        hop_ham     % Potential energy at hopping time 
        hop_frc_n   % New force at hopping time 
        hop_delta_t % time difference of hopping time and current time step 
        hop_interval% (Half-length of) interval where the local minimum is located after GD
        
        hop_pot_mat % Diabatic potentials at hopping time  
        hop_U       % Adiabatic states (eigenvectors) at hopping time
        hop_D       % Adiabatic potentials (eigenvalues) at hopping time
        hop_frc_mat % Force matrices at hopping time
        
        pos_old     % positions of particles (cell vector) of the previous step
        pos_oldold  % positions of particles (cell vector) of the previous previous step
        
        mom_old     % momenta of particles (cell vector) of the previous step
        mom_oldold  % momenta of particles (cell vector) of the previous previous step
        
        frc_old     % forces acting on "particles" (cell vector) of the previous step
        frc_oldold  % forces acting on "particles" (cell vector) of the previous previous step
        
        pot_mat_old % diabatic potential matrices of the previous step
        pot_mat_oldold
        
        frc_mat_old    % diabatic forces of the previous step
        frc_mat_oldold % diabatic forces of the previous step
        
        U_oldold       % Adiabatic states (eigenvectors) of the previous step
        
        max_round_gd  % Maximal number of rounds of gradient descent method (for lz_gr)
        acc_gd        % Accurary of gradient descent method

    end
    
    methods (Access = public)
                
        % Constructor: Set default values for class properties
        function obj = sssh_gd(lzv,n,seed)
            obj = obj@sht.sssh(lzv,n,seed);
            obj.string0 = [obj.string0 '_gd'];     % Setting the name of this class
            obj.string5 = ['With gradient descent '];
        end
        
        % Initialization
        function init_obj (obj)
            global hamilt
            
            % Inherit from initialization of superclass
            init_obj@sht.sssh ( obj );
            
            % Initialization of variables needed for GD variants
            init_gd(obj)
            
            % SSSH_GD not available in diabatic representation
            if strcmpi(hamilt.coupling.represent,'dia')
                prt.error ('SSSH_GD not available in diabatic representation')
            end
        end
        
        
        function save_previous(obj)
            
            obj.pos_oldold  = obj.pos_old;
            obj.pos_old     = obj.pos;
            
            obj.mom_oldold  = obj.mom_old;
            obj.mom_old     = obj.mom;
            
            obj.frc_oldold  = obj.frc_old;
            obj.frc_old     = obj.frc;
            
            obj.pot_mat_oldold = obj.pot_mat_old;
            obj.pot_mat_old    = obj.pot_mat;
            
            obj.frc_mat_oldold = obj.frc_mat_old;
            obj.frc_mat_old    = obj.frc_mat;
            
            obj.D_oldold  = obj.D_old;
            obj.D_old     = obj.D_new;
            
            obj.U_oldold  = obj.U_old;
            obj.U_old     = obj.U_new;
        end
        
        
        % Get probabilities of hopping from state "m" to state "n"
        function probable = prob_hop (obj,m,n,ind_m)
            
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
                
                % Call optimization (dichotomy) method
                if (~isempty(ind_c))
                    probable(ind_c) = dich_gd (obj , ind_m , ind_c , m , n);
                end
            end
        end
        
        % see separate files for the following public methods
        init_gd   ( obj )
        after_hop ( obj, mom_new,ind_h_allowed,allowed,m,n)
        probable                        = dich_gd (obj , ind_m , ind_c , m , n )
        [pos,mom,pot_mats,frc_mats,U,D] = get_quantities (obj, ind,m,n)
        
    end
        
end
    
