%--------------------------------------------------------------------------
%
% Propagate objects of class traj (trajectory bundles)
% by one substep of size time.steps.s_delta
%
%
% Velocity Verlet (third order in coordinates, first order in momenta)
% --------------------------------------------------------------------
%                        
%    q(t+dt) = q(t) + dt * p(t)/m + dt^2 * F(q(t)) / (2*m) + O(dt^4)
% 
%    p(t+dt) = p(t) + dt * ( F(q(t)) + F(q(t+dt)) ) / 2 + O(dt^2)
%
% Note that the local error order of the coordinates is 3 (not only 2)
% because of the mid point rule used in updating the momenta.
%
%    W. C. Swope, H. C. Andersen, P. H. Berens, and K. R. Wilson
%    Journal of Chemical Physics, 76(1):637–649, 1982.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019 - ... Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% see the README file for license details.

classdef v_verlet < tmp.traj.propagators.generic & handle
    
    properties (Access = private)
        frc_old     % forces at previous time step
    end
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = v_verlet
            obj.string = 'Propagator: Velocity-Verlet';
        end
        
        % Perform QC propagation: Verlet, one time step forward
        function propa (obj, state, eq_motion)
            global space time
            
            % Propagate positions by one time step
            for d = 1:space.n_dim
                state.pos{d} = state.pos{d} + time.steps.s_delta     * state.mom{d} / space.dof{d}.mass;
                state.pos{d} = state.pos{d} + time.steps.s_delta^2/2 * state.frc{d} / space.dof{d}.mass;
            end
            
            % Get (potentials and) forces at new positions
            obj.frc_old = state.frc;
            eval_V_F (eq_motion, state, time.steps.s_delta, 0)
            
            % Propagate positions by one time step
            for d = 1:space.n_dim
                state.mom{d} = state.mom{d} + time.steps.s_delta/2 * (state.frc{d}+obj.frc_old{d});
            end
                
             
        end

    end
end

