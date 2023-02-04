%--------------------------------------------------------------------------
%
% Propagate objects of class wave (wavefunctions on grids) 
% by one substep of size time.steps.s_delta
%
%
% Second order differencing
% -------------------------------------------------
%
%                               i      ^                3
%    psi(t+dt) = psi(t-dt) - 2 ---- dt H psi(t) + O ( dt ) 
%                              hbar              
%
%    A. Askar, A. S. Cakmak, J. Chem. Phys. 68, 2794 (1978)
%    http://dx.doi.org/10.1063/1.436072
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% see the README file for license details.

classdef diff_2 < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        % For future versions: add these private properties
        % psi_new
        % psi_old
    end
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = diff_2
           
        end
        
        % Display propagator, overloading default disp method
        function disp(~)
            
            prt.disp ('Finite differences, symmetrized in time           ')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp ( 'Second order differencing ' )
                        
        end

        
        % Initialize propagator
        function init (~, psi)
            
            % First step only: Initialize "old" wavefunction; no propagation yet
            psi.old = psi.dvr;
            
        end


        % Perform propagation
        function propa (~, psi,~)
            global time hamilt
            
                % Perform propagation for one substep
                if isfield(time,'pulse')
                    apply_ham(psi,[ ...
                        time.efield.grid{1}(k+time.steps.offset) ...
                        time.efield.grid{2}(k+time.steps.offset)], 0);
                else
                    apply_ham(psi,[0 0], 0);
                end
                
                for m = 1:hamilt.coupling.n_eqs
                    psi.new{m} = psi.old{m} - 2 * 1i * time.steps.s_delta * psi.new{m};
                end
                
                % Get ready for next time step
                psi.old = psi.dvr;
                psi.dvr = psi.new;
                
        end
    end
end

