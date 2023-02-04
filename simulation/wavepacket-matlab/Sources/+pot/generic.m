%--------------------------------------------------------------------------
%
% Generic properties of all potential energy class definitions
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
        dia         % Backup copy for dia=>adi transformation
        
    end
    
    methods (Access = public)
        
        % Initialize: dummy method
        function init(~)
        end
        
        % Partial copy of an instance of this class, at least a few properties
        % This is used in WavePacket in +tmp/@efield/floquet.m
        function obj2 = copy(obj1)
            obj2 = pot.generic(); % constructor of THIS class
            obj2.row = obj1.row;
            obj2.col = obj1.col;
            obj2.dvr = obj1.dvr;
        end
        
        % Display potential, overloading default disp method
        function disp (obj)
            global hamilt
            prt.disp('***************************************************************')
            prt.disp(['Potential energy for channels (' hamilt.coupling.labels{obj.row} ', ' hamilt.coupling.labels{obj.col} '):'])
        end
        
        % Grid representation of potential energy function
        function grid (obj)
            global space
            
            if ~isa (obj, 'pot.empty')
                obj.dvr = V ( obj, space.dvr );
            end
            
        end
        
    end
    
end

