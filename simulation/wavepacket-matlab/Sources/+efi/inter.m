%--------------------------------------------------------------------------
% Specific description of interpolated electric field pulse
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-... Burkhard Schmidt
%
% see the README file for license details.

classdef inter < efi.generic & handle
    
    properties (Access = public)
       file         % name of file to be opened
       method       % interpolation method
       t_conv       % Conversion factor for the time values
       f_conv       % Conversion factor for the field values
    end
    
    properties (Access = private)
       times        % tabulated time steps
       fields       % tabulated field strengths
    end
    
    methods (Access = public)
        
       % Constructor: from superclass
        function obj = inter
            obj = obj@efi.generic;
            obj.t_conv = 1;
            obj.f_conv = 1;
            obj.method = 'spline';
        end
        
       % Initialization: Read tabulate pulse
       function init ( obj )
            global time
            fid = fopen(obj.file);
            if fid == -1
                prt.error('Could not open data file')
            else
                fclose(fid);
                data = load(obj.file);
            end
            
            obj.times = data(:, 1);
            
            % It is surprisingly difficult to make matlab load complex numbers
            % from ascii files, so we work around this by loading the real and
            % imaginary values separately
            if time.efield.complex
                obj.fields = data(:,2) + 1i*data(:,3);
            else
                obj.fields = data(:, 2);
            end
            
            if obj.t_conv~=1
                obj.times = obj.times / obj.t_conv(ii);
            end
            
            if obj.f_conv~=1
                obj.fields = obj.fields / obj.f_conv(ii);
            end 
            
            obj.length = 2 * max( ...
                abs(max(obj.times(:,1))), ...
                abs(min(obj.times(:,1))));
            
            init@efi.generic ( obj );
            
        end
        
        % Display object (overwrite Matlab's default disp function)                        
        function disp(obj)
            prt.disp ('Interpolated envelope')
            prt.disp ('***************************************************************')
            prt.disp (' ')
                prt.disp(['Tabulated data from file     : '         obj.file])
                prt.disp(['Interpolation method         : '         obj.method])
                prt.disp(['Conversion factor for times  : ' num2str(obj.t_conv)])
                prt.disp(['Conversion factor for fields : ' num2str(obj.f_conv)])
            disp@efi.generic ( obj )
        end
        
        % Return interpolated of pulse
        function envelope = envelope (obj,timesteps)
            tdelayed = timesteps - obj.delay;
            envelope = interp1(obj.times, obj.fields, ...
                            tdelayed, obj.method, 0);
        end
 
        
    end
    
end