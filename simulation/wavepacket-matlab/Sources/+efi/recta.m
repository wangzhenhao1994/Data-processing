%--------------------------------------------------------------------------
% Specific description of rectangular-shaped electric field pulse
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-... Burkhard Schmidt
%
% see the README file for license details.

classdef recta < efi.generic & handle
        
    methods (Access = public)
        
       % Constructor: from superclass
        function obj = recta
            obj = obj@efi.generic;
        end
        
        % Initialization: Length of pulse (from FWHM)
        function init ( obj )
            obj.length = obj.fwhm;
            init@efi.generic ( obj );
        end

        % Display object (overload Matlab's default disp function)                
        function disp(obj)
            prt.disp ('Rectangular-shaped envelope')
            prt.disp ('***************************************************************')
            disp@efi.generic ( obj )
        end
        
        % Return envelope (shape) of pulse
        function envelope = envelope (obj,timesteps)
            tdelayed = timesteps - obj.delay;
            envelope = ones(size(tdelayed));
        end
 
        
    end
    
end