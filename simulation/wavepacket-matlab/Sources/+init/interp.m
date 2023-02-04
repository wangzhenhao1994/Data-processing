%--------------------------------------------------------------------------
% 
% Read tablulated values of wavefunctions from formatted data 
% file(s) and perform (one- or multi-dimensional) interpolation 
% Reading real and imaginary parts from two separate files
% (both of which may be omitted in which case zeros are used)
% 
% Entries for the different states:
% wav_m.dat, wav_m_im.dat for the m-th state
%
% Missing wavefunctions are ignored.
%
% Large parts of the code are taken from pot_interp.m
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% see the README file for license details.

classdef interp < handle
    
    properties (Access = public)
        pos_conv    % Conversion factor for energies
        method      % Interpolation method
        n_pts       % Number of points along each dimension
    end
    
    properties (Access = private)
        file_re     % Names of data files
        file_im     % Names of data files
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = interp
            obj.pos_conv = 1;
            obj.method = 'spline';
        end
        
        % Initialize interpolation: Set/check parameters
        function init (obj)
            global hamilt space
            
            % Check/set number of data points for multi-dimensional grids
            if space.n_dim>1
                if length(obj.n_pts)~=space.n_dim
                    error ('Inconsistent number of dimensions of tabulated data')
                end
            end
            
            % Name of input data file
            obj.file_re = cell (hamilt.coupling.n_eqs,1);
            obj.file_im = cell (hamilt.coupling.n_eqs,1);
            for m=1:hamilt.coupling.n_eqs
                obj.file_re{m} = strcat('wav_', int2str(m), '.dat');
                obj.file_im{m} = strcat('wav_', int2str(m), '_im.dat');
            end
            
        end
        
        % Display Gaussian, overloading default disp method
        function disp(obj)
            global space
            prt.disp ('Interpolation of initial wave function from tabulated data')
            prt.disp ('***********************************************************************')
            prt.disp (' ')
            prt.disp (['Conversion factor for coordinates : ' num2str(obj.pos_conv)])
            prt.disp (['Interpolation method              : '         obj.method])
             if space.n_dim>1
                prt.disp (['Number of data points in each dimension  : ' int2str(obj.n_pts)])
            end
       end
        
        % Evaluate wave function on a grid
        function wave (obj, psi)
            
            global hamilt space
            
            % Loop over all coupled states
            for m=1:hamilt.coupling.n_eqs
                if hamilt.coupling.n_eqs>1
                    prt.disp (['State: ' hamilt.coupling.labels{m}])
                end
                
                % Initialize wavefunction with zeros
                psi.dvr{m} = zeros(size(space.dvr{1}));
                 
                % Check availability of input data file (real part!)
                fid = fopen(obj.file_re{m});
                if fid == -1
                    prt.disp (['Input data file not available : ' obj.file_re{m}])
                    prt.disp ('Assuming real part of wave function to be zero');
                    prt.disp (' ')
                else
                    fclose(fid);
                    psi.dvr{m} = math.interp (...
                        space.dvr, ...
                        obj.file_re{m}, ...
                        obj.n_pts, ...
                        obj.method, ...
                        obj.pos_conv, ...
                        1);
                end
                
                % Check availability of input data file (imaginary part!)
                fid = fopen(obj.file_im{m});
                if fid == -1
                    prt.disp (['Input data file not available : ' obj.file_im{m}])
                    prt.disp ('Assuming imaginary part of wave function to be zero');
                    prt.disp (' ')
                else
                    fclose(fid);
                    psi.dvr{m} = psi.dvr{m} + ...
                        1i * math.interp (...
                        space.dvr, ...
                        obj.file_im{m}, ...
                        obj.n_pts, ...
                        obj.method, ...
                        obj.pos_conv, ...
                        1);
                end
                
                prt.disp(' ')
            end

        end
        
    end
end