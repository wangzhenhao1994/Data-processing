%--------------------------------------------------------------------------
%
% Generic properties of all negative imaginary potential class definitions
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.


classdef generic < handle
    
    properties (Access = public) 
        
        ind         % Index of (coupled) channel
        
        dvr         % Grid representation (in N dimensions)

    end

    methods (Access = public)
        
        % Initialize: dummy method
        function init(~)
        end
        
        % Display potential, overloading default disp method
        function disp (obj)
            prt.disp('***************************************************************')
            prt.disp(['Negative imaginary potential for channel (' int2str(obj.ind) '):'])
        end
        
        % Grid representation of negative imaginary potential
        function grid (obj)
            global space
            
           if ~isa (obj, 'nip.empty')
                obj.dvr = W ( obj, space.dvr );
           end
            
        end
        
    end

end

