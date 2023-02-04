%--------------------------------------------------------------------------
%
% Generic properties of all polarizability class definitions
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
        p_1         % polarization direction
        p_2         % polarization direction
        
        row         % row index (in diabatic representation)
        col         % column index (in diabatic representation)
        
        dvr         % Grid representation (in N dimensions)

    end
    
    methods (Access = public)
        
        % Initialize: dummy method
        function init(~)
        end
        
        % Display polarizability, overloading default disp method
        function disp (obj)
            global hamilt
            prt.disp('***************************************************************')
            prt.disp([ int2str(obj.p_1) '-' int2str(obj.p_2) '-th component of polarizability for channels (' hamilt.coupling.labels{obj.row} ', ' hamilt.coupling.labels{obj.col} '):'])
        end

        % Grid representation of polarizability
        function grid (obj)
            global space
            
            if ~isa (obj, 'pol.empty')
                obj.dvr = alpha ( obj, space.dvr );
            end
            
        end
        
    end
    
end

