%--------------------------------------------------------------------------
%
% Generic properties of all system-bath coupling class definitions
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
        
        row         % row index (in diabatic potential matrix)
        col         % column index (in diabatic potential matrix)
        
        dvr         % Grid representation (in N dimensions)

    end
    
    methods (Access = public)

        % Initialize: dummy method
        function init(~)
        end
        
        % Display system-bath coupling, overloading default disp method
        function disp (obj)
            global hamilt
            prt.disp('***************************************************************')
            prt.disp(['System-bath coupling for channels (' hamilt.coupling.labels{obj.row} ', ' hamilt.coupling.labels{obj.col} '):'])
        end

        % Grid representation of system-bath coupling
        function grid (obj)
            global space
            
            if ~isa (obj, 'sbc.empty')
                obj.dvr = chi ( obj, space.dvr );
            end
            
        end
        
    end
    
end

