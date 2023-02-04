%--------------------------------------------------------------------------
%
% "Single switch surface hopping" class definitions
% using the GD="gradient descent" method to exactly 
% MC variant: ... text to be added by Leo ...
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-.... Leonardo Cancissu Araujo
%
% see the README file for license details.

classdef sssh_gd_mc < sht.sssh_gd & handle
    
    properties (Access = public)
        hop_direction
    end
    
    methods (Access = public)
        
        % Constructor: Set default values for class properties
        function obj = sssh_gd_mc(lzv,n,seed)
            obj = obj@sht.sssh_gd(lzv,n,seed);
            obj.string0 = [obj.string0 '_mc'];  
            obj.string5 = 'With gradient descent and with multiple channels';
        end

        
        % Initialization
        function init_obj (obj)
            global hamilt
            
            % Inherit from initialization of superclass
            init_obj@sht.sssh_gd ( obj );
            
            % LZGD Multiple Channels not available in diabatic representation
            if strcmpi(hamilt.coupling.represent,'dia')
                prt.error ('LZGD Multiple Channels is not available in diabatic representation')
            end
            
            obj.hop_direction             = zeros(obj.n_p,1);
        end
                
        
        function ind_hop = get_ind_hop(obj , ind_m , m , n)
            
            % Get probabilities of hopping (from sub-classes)
            probable = prob_hop (obj,m,n,ind_m);
            
            % If probabilities are negative, set them to zero
            probable ( probable<0 ) = 0;
            
            % Uniform random numbers in the interval (0,1)
            zeta = rand(size(ind_m));
            
            % Find indices of hopping trajectories by comparison with
            % zeta, a uniform random number in the interval (0,1)
            ind_hop = ind_m (zeta>=0 & zeta<probable );
        end
        
        
        
        % see separate files for the following public methods
        traj_hop  ( obj )                % Perform the trajectory hopping
        after_hop ( obj, mom_new,ind_h_allowed,allowed,m,n)
        [ind_hop_allowed,ind_hop_forbidden] = perform_hop (obj, ind_hop, m, n)
        [pos, mom, frc, pot_mats, frc_mats, D_1, U_1] = verlet_prop_back(obj,pos,mom,frc,n)
    end

end
        
    
