%--------------------------------------------------------------------------
%
% Fully classical trajectories
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

classdef classical < handle
    
    properties (Access = public)
    end    

    methods (Access = public)
        
        % Construct object
        function obj = classical
        end
            
        % Initialize
        function eq_init (obj, state, propa_type)   
            init ( propa_type, state, obj)
        end
        
        % Propagate
        function eq_propa (obj, state, propa_type) 
            propa ( propa_type, state, obj)
        end
        
        % Optionally solve TDSEs attached to trajectories
        function solve_tdse(obj, state, step_size, D, U, U_old)
            if (state.extra_tdse_solving)
                tdse_adi (state, step_size, D, U, U_old)
            end
        end
        
        % Calculate forces in adiabatic representation
        function frc_adi = eq_frc (obj, pot_mats, frc_mats, D, U, m, psi) 
            frc_adi = ham.frc_adi (pot_mats,frc_mats,U,m);
        end
        
        % see separate files for the following public methods
        eval_V_F  ( obj, state, step_size, init)% Evaluate potential energy and forces 
        
    end
end



    
