%-------------------------------------------------------------------------
%
% Propagate objects of class traj (trajectory bundles) 
% by one substep of size time.steps.s_delta
%
% Leap frog (second order):
% -------------------------
%
%    p(t+dt/2) = p(t     ) + dt * F(q(t)    )/2 
%    q(t+dt  ) = q(t     ) + dt *   p(t+dt/2)/m
%    p(t+dt  ) = p(t+dt/2) + dt * F(q(t+dt  ))/2 
%
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-20.. Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

classdef leap_frog < tmp.traj.propagators.generic & handle   

    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = leap_frog
            obj.string = 'Propagator: Leap frog (Strang-Marchuk)';
        end
        
        function propa(~,state,eq_motion)
            global space time
                
                % Propagate momenta by half time step
                for d = 1:space.n_dim
                    state.mom{d} = state.mom{d} + 0.5 * time.steps.s_delta * state.frc{d};
                end
                
                % Propagate positions by full time step
                for d = 1:space.n_dim
                    state.pos{d} = state.pos{d} + time.steps.s_delta * state.mom{d} / space.dof{d}.mass;
                end
                
                % Get potentials and forces at new positions
                eval_V_F (eq_motion, state, time.steps.s_delta, 0)
                
                % Propagate momenta by half time step
                for d = 1:space.n_dim
                    state.mom{d} = state.mom{d} + 0.5 * time.steps.s_delta * state.frc{d};
                end
        end 
        
    end
end


