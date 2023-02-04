%--------------------------------------------------------------------------
% Specific description of Gaussian-shaped electric field pulse
%
%                    (        ( time - delay )^2      )
%     f (time) = exp ( -  --------------------------- )
%                    (    2 ( fwhm / sqrt(8*ln2) )^2  )
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-... Burkhard Schmidt
%
% see the README file for license details.

classdef gauss < efi.generic & handle
    
    methods (Access = public)
        
        % Constructor: from superclass
        function obj = gauss
            obj = obj@efi.generic;
        end
        
        % Initialization: Length of pulse (from FWHM)
        function obj = init ( obj )
            obj.length = obj.fwhm * sqrt(log(eps)/log(1/2));
            init@efi.generic ( obj );
        end
        
        % Display object (overload Matlab's default disp function)
        function disp(obj)
            prt.disp ('Gaussian-shaped envelope')
            prt.disp ('***************************************************************')
            disp@efi.generic ( obj )
        end
        
        % Return envelope (shape) of pulse
        function envelope = envelope (obj,timesteps)
            tdelayed = timesteps - obj.delay;
            envelope = exp ( - tdelayed.^2 / ...
                ( 2 * (obj.fwhm/sqrt(8*log(2)))^2) );
        end
 
        
    end
    
end