%--------------------------------------------------------------------------
% 
% Read tabulated values of dipole moment functions from formatted data 
% file(s) and perform (one- or multi-dimensional) interpolation 
% 
% Entries of dipole moment matrices:
% Permanent dipole moments: 'd_x_m.dat', 'd_y_m.dat' with m=1,2,...
% Transition dipole moments: 'd_x_m_n.dat', 'd_y_m_n.dat' with n>m
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef interp < dip.generic & handle
    
    properties (Access = public)
        
        pos_conv    % Conversion factor for coordinates
        dip_conv    % Conversion factor for energies
        method      % Interpolation method
        n_pts       % Number of points along each dimension
        
        file        % Nsme of data file

    end
    
    methods (Access = public)
        
       % Constructor: Set default values
        function obj = interp
            obj.pos_conv = 1;
            obj.dip_conv = 1;
            obj.method = 'spline';
        end
        
        % Initialize dipole moment: Set/check parameters
        function init (obj)
            
            % Check/set number of data points for multi-dimensional grids
            global space
            if space.n_dim>1
                if length(obj.n_pts)~=space.n_dim
                    prt.error ('Inconsistent number of dimensions of tabulated data')
                end
            end
            
            % Polarization
            switch obj.pol
                case 1
                    p = 'x';
                case 2
                    p = 'y';
                case 3
                    p = 'z';
                otherwise
                    prt.error ('Polarization direction has to be 1,2,3 for x,y,z')
            end
            
            % Name of input data file
            if obj.row==obj.col
                obj.file = strcat('d_', p, '_', int2str(obj.row), '.dat');
            else
                obj.file = strcat('d_', p, '_', int2str(obj.row), '_', int2str(obj.col), '.dat');
            end
            
            % Check availability of input data file
            fid = fopen(obj.file);
            if fid == -1
                    prt.error (['Input data file missing : ' obj.file])
            else
                fclose(fid);
            end

        end
        
        % Display dipole moment, overloading default disp method
        function disp (obj)
            global hamilt space
            disp@dip.generic(obj)
            prt.disp ('Interpolating tabulated values of dipole moment')
            prt.disp ('***************************************************************')
            prt.disp ( ' ' )
            prt.disp (['Conversion factor for coordinates : ' num2str(obj.pos_conv)])
            prt.disp (['Conversion factor for dipoles     : ' num2str(obj.dip_conv)])
            prt.disp (['Interpolation method              : '         obj.method])
            if space.n_dim>1
                prt.disp (['Number of data points in each dimension  : ' int2str(obj.n_pts)])
            end
            if obj.row==obj.col
                prt.disp (['Read dipole moment (' hamilt.coupling.labels{obj.row} ') from data file : ' obj.file])
            else
                prt.disp (['Read dipole moment (' hamilt.coupling.labels{obj.row},', ',hamilt.coupling.labels{obj.col} ') from data file : ' obj.file])
            end
            prt.disp (' ')
        end
        
        % Evaluate dipole moment function
        function mu = mu(obj,r)
             
            mu = math.interp (...
                    r, ...
                    obj.file, ...
                    obj.n_pts, ...
                    obj.method, ...
                    obj.pos_conv, ...
                    obj.dip_conv);

        end
    end
end