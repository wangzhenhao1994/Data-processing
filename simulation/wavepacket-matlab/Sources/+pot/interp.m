%--------------------------------------------------------------------------
% 
% Read tabulated values of potential energy functions from formatted data 
% file(s) and perform (one- or multi-dimensional) interpolation 
% 
% Entries of diabatic potential energy matrices:
% Diagonal elements: 'pot_m.dat' with m=1,2,...
% Off-diagonal elements: 'pot_m_n.dat' with n>m
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef interp < pot.generic & handle
    
    properties (Access = public)
        
        pos_conv    % Conversion factor for coordinates
        pot_conv    % Conversion factor for energies
        method      % Interpolation method
        n_pts       % Number of points along each dimension
        
        file        % Name of data file
        
    end
    
    methods (Access = public)
        
       % Constructor: Set default values
        function obj = interp
            obj.pos_conv = 1;
            obj.pot_conv = 1;
            obj.method = 'spline';
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            global space
            
            % Check/set number of data points for multi-dimensional grids
            if space.n_dim>1
                if isempty(obj.n_pts)
                    prt.error ('Number of tabulated data has to be specified')
                end
                if length(obj.n_pts)~=space.n_dim
                    prt.error ('Inconsistent number of dimensions of tabulated data')
                end
            end
            
            % Name of input data file
            if obj.row==obj.col
                obj.file = strcat('pot', '_', int2str(obj.row), '.dat');
            else
                obj.file = strcat('pot', '_', int2str(obj.row), '_', int2str(obj.col), '.dat');
            end
            
            % Check availability of input data file
            fid = fopen(obj.file);
            if fid == -1
                    prt.error (['Input data file missing : ' obj.file])
            else
                fclose(fid);
            end

        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            global space
            disp@pot.generic (obj)
            prt.disp ('Interpolating tabulated values of potential energy')
            prt.disp ('***************************************************************')
            prt.disp ( ' ' )
            prt.disp (['Conversion factor for coordinates   : ' num2str(obj.pos_conv)])
            prt.disp (['Conversion factor for energies      : ' num2str(obj.pot_conv)])
            prt.disp (['Interpolation method                : '         obj.method])
            prt.disp (['Read tabulated potential from file  : '         obj.file])
            if space.n_dim>1
                prt.disp (['Number of points in each dimension  : ' int2str(obj.n_pts)])
            end
            prt.disp (' ')
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)
            
            V = math.interp (...
                    r, ...
                    obj.file, ...
                    obj.n_pts, ...
                    obj.method, ...
                    obj.pos_conv, ...
                    obj.pot_conv);

        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            prt.error ('Code for calculation of forces still missing')
        end
        
    end
end