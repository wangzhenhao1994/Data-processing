%--------------------------------------------------------------------------
% Specific description of Sine-squared-shaped electric field pulse
% 
%                    2 (     time - delay    pi )
%     f  (time) = sin  ( pi ------------- + --- )
%                      (        2 fwhm       2  )
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-... Burkhard Schmidt
%
% see the README file for license details.

classdef sin_2 < efi.generic & handle
    
    methods (Access = public)
        
       % Constructor: from superclass
        function obj = sin_2
            obj = obj@efi.generic;
        end
        
        % Initialization: Length of pulse (from FWHM)
        function init ( obj )
            obj.length = 2* obj.fwhm;
            init@efi.generic ( obj );
        end

        % Display object (overload Matlab's default disp function)     
        function disp(obj)
            prt.disp ('Sine-squared-shaped envelope')
            prt.disp ('***************************************************************')
            disp@efi.generic ( obj )
        end
        
        % Return envelope (shape) of pulse
       function envelope = envelope (obj,timesteps)
            tdelayed = timesteps - obj.delay;
            envelope = cos ( pi/2 * tdelayed/obj.fwhm ).^2;
        end
 
        
    end
    
end