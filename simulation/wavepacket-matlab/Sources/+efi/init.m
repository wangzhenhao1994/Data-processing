%--------------------------------------------------------------------------
%
% The external electric field is assumed to be a sequence of (overlapping 
% or non-overlapping) pulses, thereby mimicking modern time-resolved 
% spectrosopic experiments employing (ultra-)short laser pulses.  
%
% Each pulse has a constant carrier frequency modulated by an envelope of 
% different shape. Furthermore it is characterized by the field amplitude,
% polarization, time delay (center of pulse), duration(FWHM), and phase
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function init
global time

%% Set the number of pulses to default values and start output
if isfield(time, 'pulse')

    if ~isfield(time,'efield') % If not already existing
        time.efield = tmp.efield; % Construct object
    end
    
    prt.disp ('***************************************************************')
    prt.disp ('External electric field as a sequence of pulses')
    prt.disp ('***************************************************************')
    prt.disp ( ' ' )
    prt.disp ( [ 'Using (Floquet) dressed states  : ' int2str(time.efield.dressed) ] )
    prt.disp ( [ 'Using complex-valued fields     : ' int2str(time.efield.complex) ] )
    prt.disp ( [ 'Number of pulses                : ' int2str(length(time.pulse) ) ] )
    prt.disp ( ' ' )
        
    %% Set up individual pulses
    for p=1:length(time.pulse)
        
        prt.disp ('***************************************************************')
        prt.disp ( [ 'Parameters of electric field pulse : ' int2str(p) ] )
        
        if ~isempty (time.pulse{p})
            time.pulse{p}.ind = p; % Tell each pulse its index
            init ( time.pulse{p} ); % Initialize pulse
            disp ( time.pulse{p} ); prt.disp ( ' ' ) % Display pulse
        else
            prt.disp ( 'Not available' )
            prt.disp ('***************************************************************')
            prt.disp ( ' ' )
            
        end
    end
    
    % Initialize fields
    init (time.efield);
    
    % Optionally switch to Floquet dressed
    if time.efield.dressed
        floquet (time.efield);
    end
    
    
else
    
    prt.disp ('***************************************************************')
    prt.disp ('No electric field available                                    ')
    prt.disp ('***************************************************************')
    prt.disp (' ')
end

end




