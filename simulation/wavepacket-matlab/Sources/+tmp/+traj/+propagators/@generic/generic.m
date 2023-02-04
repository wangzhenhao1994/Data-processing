%--------------------------------------------------------------------------
%
% Generic class for all "traj propagators" class definitions
%
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

classdef generic < handle
    
    properties (Access = public)
        string      % Descriptive name of this propagator
    end    

    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = generic
        end
        
        % Display propagator, overloading default disp method
        function disp(obj)
            prt.disp ('***************************************************************')
            prt.disp ('Numeric propagation scheme:                       ')
            prt.disp (obj.string)
            prt.disp ('***************************************************************')
            prt.disp (' ')
        end
            
        % Initialize QC propagator
        function init (obj, state, eq_motion)     
            
            % Get initial energies
            eval_V_F  ( eq_motion, state, [], 1);
            apply_ham ( state );
            
            disp(obj);
            
        end
        
    end
end



    
