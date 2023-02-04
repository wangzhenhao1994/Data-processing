%-------------------------------------------------------------------------
%
% Propagate objects of class ket (state vectors) or
% or class rho (density matrices) by one time step 
% of size tau using a Runge-Kutta integrator
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-20.. Burkhard Schmidt's group
%
% see the README file for license details.

classdef RungeKutta < handle
    
    properties (Access = public)
        integr      % Handle of integrator
        order       % Error order of integrator
    end
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = RungeKutta (order)
            if isempty (order)
                obj.order = 6;
            end
        end
              
        % Display propagator, overloading default disp method
        function disp(obj)
            prt.disp ('Using Runge-Kutta ODE solvers')
            prt.disp ('***************************************************************')
            prt.disp (' ') 
            prt.disp (['Order of propagator : ',  int2str(obj.order) ]) 
        end
        
        % Initialize propagator
        function init (obj)
            switch obj.order
                case 4
                    obj.integr = 'tmp.ode.RungeKutta.RuKu4';
                case 5
                    obj.integr = 'tmp.ode.RungeKutta.RuKu5';
                case 6
                    obj.integr = 'tmp.ode.RungeKutta.RuKu6';
                otherwise
                    prt.error ('Wrong choice of error order')
            end
        end
        
        % Perform forward  propagation
        function forward (obj,step,derivative)
            global control
            if derivative
                control.x.forward(:,step) = feval ( ...
                    obj.integr, ...
                    @tmp.ode.rhs_x, ...
                    control.x.forward(:,step-1), ...
                    control.u.forward(:,step-1), ...
                    control.d.forward(:,step-1), ...
                   +control.t.delta);
            else
                control.x.forward(:,step) = feval ( ...
                    obj.integr, ...
                    @tmp.ode.rhs_x, ...
                    control.x.forward(:,step-1), ...
                    control.u.forward(:,step-1), ...
                    0, ...
                   +control.t.delta);
            end
        end
               
        % Perform backward propagation
        function backward (obj,step)
            global control
                control.x.backward(:,step) = feval ( ...
                    obj.integr, ...
                    @tmp.ode.rhs_z, ...
                    control.x.backward(:,step+1), ...
                    control.u.backward(:,step+1), ...
                    control.d.backward(:,step+1), ...
                   -control.t.delta);
        end
        
    end

    
    methods (Static)
        xout = RuKu4(RHS,xin,u,dudt,tau)
        xout = RuKu5(RHS,xin,u,dudt,tau)
        xout = RuKu6(RHS,xin,u,dudt,tau)
    end
end


