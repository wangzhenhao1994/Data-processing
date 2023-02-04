%--------------------------------------------------------------------------
%
% Propagate objects of class traj (trajectory bundles)
% by one substep of size time.steps.s_delta
%
% Stoermer-Verlet (third order)
% -----------------------------
%                                   F(q(t))   2        4
%    q(t+dt) = - q(t-dt) + 2 q(t) + ------- dt + O ( dt ) 
%                                    m
%
% Although they are not part of the propagation scheme,
% knowledge of the momenta is sometimes required, e.g.
% for calculating the kinetic energy. In those cases,
% the momenta are obtained from the following relation
%
%   p(t) = [ q(t+dt) - q(t-dt) ] / (2*dt) + O(dt^2)
%
% Note that the local error order of the coordinates is 3 (not only 2)
% because of the symmetrization in time.
% 
%    L. Verlet, Phys. Rev. 159, 98-103 (1967)  
%    http://dx.doi.org/10.1103/PhysRev.159.98
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019 - ... Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% see the README file for license details.

classdef s_verlet < tmp.traj.propagators.generic & handle
    
    properties (Access = private)
        pos_old     % "old" positions
        pos_new     % "new" positions
    end
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = s_verlet
            obj.string = 'Propagator: Stoermer-Verlet';
        end

        % Initialize QC propagator: Verlet, one time step backward
        function init (obj, state, eq_motion)
            global space time
            
            % Inherit from "generic" super class
            init@tmp.traj.propagators.generic(obj, state, eq_motion);
            
            % Propagate positions one step backward in time (first order)
            for d = 1:space.n_dim
                obj.pos_old{d} = state.pos{d} - time.steps.s_delta * state.mom{d} / space.dof{d}.mass;
            end
                        
        end
        
        % Perform QC propagation: Verlet, one time step forward
        function propa (obj, state, eq_motion)
            global space time
            
            % Get (potentials and) forces
            eval_V_F (eq_motion, state, time.steps.s_delta, 0)
            
            % Propagate positions by one time step
            % Get momenta by simple differencing
            for d = 1:space.n_dim
                obj.pos_new{d} = 2*state.pos{d} - obj.pos_old{d} + time.steps.s_delta^2 * state.frc{d} / space.dof{d}.mass;
                state.mom{d} = ( obj.pos_new{d} - obj.pos_old{d} ) / (2*time.steps.s_delta) * space.dof{d}.mass;
                obj.pos_old{d} = state.pos{d};
                state.pos{d} = obj.pos_new{d};
            end
            
        end
        
    end
end

