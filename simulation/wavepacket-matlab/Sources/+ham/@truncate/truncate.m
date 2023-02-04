%--------------------------------------------------------------------------
%
% Truncating the grid representation of the Hamiltonian operator
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef truncate < handle
    
    properties (Access=public)
        
        e_min       % truncate energies below this value
        e_max       % truncate energies above this value
        delta       % range of truncated energies
        
    end
    
    methods (Access=public)
        
        % Constructor: Set default values
        function obj = truncate
            obj.e_min = -Inf;              % truncate energies below this value
            obj.e_max = +Inf;              % truncate energies above this value
        end
        
        % Initialization
        function init (obj)
            obj.delta = obj.e_max - obj.e_min;
        end
        
        % Display properties of object
        function disp (obj)
            
            prt.disp ('***************************************************************')
            prt.disp ('Initialize truncation of Hamiltonian     ')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp ( [ 'Lower truncation of energies  : ' num2str(obj.e_min) ] )
            prt.disp ( [ 'Upper truncation of energies  : ' num2str(obj.e_max) ] )
            prt.disp ( [ 'Truncation range              : ' num2str(obj.delta) ] )
            prt.disp ( ' ' )
            
        end
        
        % Truncating potential and kinetic energy
        function trunc_pot_kin (obj)
            global hamilt space
            
            % Manually truncate potential energy
            for m = 1:hamilt.coupling.n_eqs
                if ~isempty (hamilt.pot{m,m}.dvr)
                    hamilt.pot{m,m}.dvr(hamilt.pot{m,m}.dvr>obj.e_max) = obj.e_max;
                    hamilt.pot{m,m}.dvr(hamilt.pot{m,m}.dvr<obj.e_min) = obj.e_min;
                end
            end
            
            % Estimate spectral range of Hamiltonian (neglect off-diagonal coupling)
            hamilt.kin_max = 0;
            for k = 1:space.n_dim
                hamilt.kin_max = hamilt.kin_max + min(space.dof{k}.kin_max, obj.delta);
            end
            
            % Estimate spectral range of additional kinetic operators
            if isfield(hamilt, 'kin')
                for k = 1:length(hamilt.kin)
                    hamilt.kin_max = hamilt.kin_max + min(hamilt.kin{k}.kin_max, obj.delta);
                end
            end
            
            % Estimate spectral range of potential energy
            hamilt.pot_min = +realmax; % Largest positive floating-point number
            hamilt.pot_max = -realmax;
            for m=1:hamilt.coupling.n_eqs
                if ~isempty (hamilt.pot{m,m}.dvr)
                    hamilt.pot_min = min ( hamilt.pot_min, min ( hamilt.pot{m,m}.dvr(:) ) );
                    hamilt.pot_max = max ( hamilt.pot_max, max ( hamilt.pot{m,m}.dvr(:) ) );
                else
                    hamilt.pot_min = min ( hamilt.pot_min, 0 );
                    hamilt.pot_max = max ( hamilt.pot_max, 0 );
                end
            end
            
            % Set spectral range of Hamiltonian manually or from min/max of kin/pot
            if ~isfield(hamilt, 'range')
                hamilt.range = [];
            end
            if ~isfield(hamilt.range, 'e_min')
                hamilt.range.e_min = hamilt.pot_min;
            end
            if ~isfield(hamilt.range, 'e_max')
                hamilt.range.e_max = hamilt.kin_max + hamilt.pot_max;
            end
            hamilt.range.delta   = hamilt.range.e_max - hamilt.range.e_min;
            
            %% Output
            prt.disp ('***************************************************************')
            prt.disp ('Spectral range of Hamiltonian     ')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp ( [ 'Maximum of kinetic energy     : ' num2str(hamilt.kin_max) ] )
            prt.disp ( [ 'Minimum of potential energy   : ' num2str(hamilt.pot_min) ] )
            prt.disp ( [ 'Maximum of potential energy   : ' num2str(hamilt.pot_max) ] )
            prt.disp ( ' ' )
            prt.disp ( [ 'Minimum of grid Hamiltonian   : ' num2str(hamilt.range.e_min) ] )
            prt.disp ( [ 'Maximum of grid Hamiltonian   : ' num2str(hamilt.range.e_max) ] )
            prt.disp ( [ 'Spectral range of Hamiltonian : ' num2str(hamilt.range.delta) ] )
            prt.disp ( ' ' )
            
            
        end
        
        
    end
    
end

