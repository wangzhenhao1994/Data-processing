%--------------------------------------------------------------------------
%
% Generic properties of all additional multiplicative operators class definitions
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
        label       % Label of this AMO function
        ind         % Index of this AMO function
        
        dvr         % Grid representation (in N dimensions)
        
    end
  
    methods (Access = public)
        
        % Initialize: dummy method
        function init(~)
        end
        
       % Display AMO, overloading default disp method
        function disp (obj)
            prt.disp('***************************************************************')
            prt.disp([ int2str(obj.ind) '-th additional multiplicative operator'])
        end
 
        % Grid representation of additional multiplicative operators
        function grid (obj)
            global space
            
            if ~isa (obj, 'amo.empty')
                obj.dvr = A ( obj, space.dvr );
            end

        end
        
    end
 
end

