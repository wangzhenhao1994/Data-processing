%--------------------------------------------------------------------------
% General description of an electric field pulse as Matlab objects
%--------------------------------------------------------------------------

classdef generic < handle
    
    properties (Access = public)
        
        fwhm        % FWHM: full width at half maximum
        length      % Full duration of pulse
        delay       % Time delay (center of pulse)
        
        polar       % Polarization angle [rad]
        ampli       % Field amplitude
        strength    % Same, but components along x, y
        
        frequ       % Frequency (photon energy)
        linear      % Linear chirp
        quadratic   % Quadratic chirp
        phase       % Phase shift
        
        ind         % index of this pulse
        
    end
    
    properties  (Access = private)
        l_cycle     % length of optical cycles
        n_cycle     % number of optical cycles
        
    end
    
    methods  (Access = public)
        
        % Constructor: Set default values
        function obj = generic
            obj.fwhm = 0;
            obj.length = 0;
            obj.delay = 0;
            obj.polar = 0;
            obj.ampli = 0;
            obj.frequ = 0;
            obj.linear = 0;
            obj.quadratic = 0;
            obj.phase = 0;
        end
        
        % Initialize some properties of the pulse
        function init ( obj )
            obj.strength(1) = obj.ampli * cos ( obj.polar );
            obj.strength(2) = obj.ampli * sin ( obj.polar );
            if obj.frequ~=0
                obj.l_cycle = 2*pi/obj.frequ;
                obj.n_cycle = obj.length / obj.l_cycle;
            end
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            prt.disp ( ' ' )
            prt.disp (['FWHM duration of pulse      : ' num2str(obj.fwhm)])
            prt.disp (['Full duration of pulse      : ' num2str(obj.length)])
            prt.disp (['Time delay (center of pulse): ' num2str(obj.delay)])
            prt.disp ( ' ' )
            prt.disp (['Polarization angle [rad]    : ' num2str(obj.polar)])
            prt.disp ( ' ' )
            prt.disp (['Field amplitude             : ' num2str(obj.ampli)])
            prt.disp (['Component along x           : ' num2str(obj.strength(1))])
            prt.disp (['Component along y           : ' num2str(obj.strength(2))])
            prt.disp ( ' ' )
            prt.disp (['Frequency (photon energy)   : ' num2str(obj.frequ)])
            prt.disp (['Linear chirp                : ' num2str(obj.linear)])
            prt.disp (['Quadratic chirp             : ' num2str(obj.quadratic)])
            prt.disp (['Phase shift                 : ' num2str(obj.phase)])
            if obj.frequ~=0
                prt.disp ( ' ' )
                prt.disp (['Duration of optical cycle   : ' num2str(obj.l_cycle)])
                prt.disp (['Number of optical cycles    : ' num2str(obj.n_cycle)])
            end
        end
        
        % Calculate the oscillatory portion of the electric field
        % Constant or (linearly or quadratically) chirped frequencies
        function oscillate = oscillate (obj,timesteps)
            tdelayed = timesteps - obj.delay;
            omega = obj.frequ * ones(size(tdelayed));
            omega = omega + obj.linear * tdelayed + obj.quadratic/2 * tdelayed.^2;
            oscillate = cos ( omega .* tdelayed + obj.phase );
        end
        
    end
    
end


