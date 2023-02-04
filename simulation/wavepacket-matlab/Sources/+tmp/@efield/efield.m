%------------------------------------------------------------------------------
%
% Temporal discretization: Dealing with electric fields
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2010-2011 Ulf Lorenz
%
% see the README file for license details.

classdef efield < handle
    
    properties (Access=public)
        
        n_pulse     % number of pulses
        dressed     % toggle deressed vs. bare states
        photons     % dressed with this many photons
        complex     % toggle use of complex-valued fields
        
        has_x       % Check x directions
        has_y       % Check y directions
        
        max_ampli   % Maximum amplitude
        
        grid        % values of el. field for all time (sub!) steps
        
    end
    
    methods (Access=public)
        
        % Constructor: Set trivial default values
        function obj = efield
            
            obj.n_pulse = 0;
            obj.complex = false;
            obj.dressed = false;
            
            obj.has_x = false;
            obj.has_y = false;
            
            obj.max_ampli = 0;
            
        end
        
        % Initialization: Check / set parameters
        function init (obj)
            
            global hamilt time
            
            % Error checking
            if time.efield.complex && hamilt.coupling.n_eqs > 2
                prt.error('Cannot use complex electric fields with more than 2 electronic states');
            end
            
            % Number of electric field pulses
            obj.n_pulse = length(time.pulse);
            
            % Max amplitude
            for p=1:length(time.pulse)
                if ~isempty(time.pulse{p})
                    if abs(time.pulse{p}.ampli) > abs(obj.max_ampli)
                        obj.max_ampli = time.pulse{p}.ampli;
                    end
                end
            end
            
            % Check x|y directions
            for p=1:length(time.pulse)
                if ~isempty(time.pulse{p})
                    if time.pulse{p}.polar~=pi/2
                        obj.has_x = true;
                    end
                    if time.pulse{p}.polar~=0
                        obj.has_y = true;
                    end
                end
            end
            
            % Set up grid for values of the electric field
            obj.grid{1}    = zeros(length(time.steps.s_grid),1);
            obj.grid{2}    = zeros(length(time.steps.s_grid),1);
            
            % Evaluate electric field for all time (sub!) steps
            obj.grid = efi.eval (time.steps.s_grid);
            
        end
        
        % Floquet function: see extra file
        floquet (obj)
        
    end
    
end
