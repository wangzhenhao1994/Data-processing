%-------------------------------------------------------------------------
%
% Propagate objects of class ket (state vectors) or
% or class rho (density matrices) by one main step 
% of size time.steps.m_delta using one of Matlab's 
% builtin ODE solvers. Internally these ODE solvers 
% may use finer grids, typically chosen adaptively
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-20.. Burkhard Schmidt's group
%
% see the README file for license details.

classdef matlab < handle
    
    properties (Access = public)
        integr      % Handle of integrator
        reltol      % Relative tolerance
        options     % Options structure for Matlab's ODE solvers
    end
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = matlab (integr,reltol)
            if isempty (integr)
                obj.integr = 'ode45';
            end
            if isempty (reltol)
                obj.reltol = 1e-6;
            end
        end
              
        % Display propagator, overloading default disp method
        function disp(obj)
            prt.disp ('Using MATLAB''s built-in adaptive ODE solvers')
            prt.disp ('***************************************************************')
            prt.disp (' ') 
            prt.disp (['Choice of propagator : ',         obj.integr ]) 
            prt.disp (['Relative tolerance   : ', num2str(obj.reltol)]) 
        end
        
        % Initialize propagator
        function init (obj, ~)            
            obj.options = odeset ('reltol',obj.reltol);
        end

        % Perform propagation
        function propa (obj, state, step)
            global time
            
            [~,x_new] = feval ( ...
                obj.integr, ...
                @tmp.ode.evolve, ...
                [time.steps.m_grid(step-1) time.steps.m_grid(step)], ...
                state.x, ...
                obj.options );
            state.x = x_new(end,:).';
            
        end
        
    end
end


