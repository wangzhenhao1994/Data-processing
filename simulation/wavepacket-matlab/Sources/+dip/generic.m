%--------------------------------------------------------------------------
%
% Generic properties of all dipole moment class definitions
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
        pol         % polarization direction
        
        row         % row index (in diabatic representation)
        col         % column index (in diabatic representation)
        
        dvr         % Grid representation (in N dimensions)
        
    end
    
    methods (Access = public)
        
        % Initialize: dummy method
        function init(~)
        end
        
        % Partial copy of an instance of this class, at least a few properties
        % This is used in WavePacket in +tmp/@efield/floquet.m only
        function obj2 = copy(obj1)
            obj2 = dip.generic(); % constructor of THIS class
            obj2.pol = obj1.pol;
            obj2.row = obj1.row;
            obj2.col = obj1.col;
            obj2.dvr = obj1.dvr;
        end
        
        % Display dipole moment, overloading default disp method
        function disp (obj)
            global hamilt
            prt.disp('***************************************************************')
            prt.disp([ int2str(obj.pol) '-th component of dipole moment for channels (' hamilt.coupling.labels{obj.row} ', ' hamilt.coupling.labels{obj.col} '):'])
        end
        
        % Grid representation of dipole moment
        function grid (obj)
            global space
            
            if ~isa (obj, 'dip.empty')
                obj.dvr = mu ( obj, space.dvr );
            end
            
        end
        
    end
    
end

