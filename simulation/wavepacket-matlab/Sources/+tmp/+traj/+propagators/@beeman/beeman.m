%-------------------------------------------------------------------------
%
% Propagate objects of class traj (trajectory bundles) 
% by one substep of size time.steps.s_delta
%
% Beeman (third order):
% ---------------------
%
%     Reference: D. Beeman
%     J. Comp. Phys. 20, 130 (1976)
%     DOI:10.1016/0021-9991(76)90059-0
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

classdef beeman < tmp.traj.propagators.generic & handle
    
    properties (Access = private)
        frc_old     % Forces from previous time step
    end    
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = beeman
            obj.string = 'Propagator : Beeman (3rd order)';
            obj.frc_old = [];
        end

        % Perform QC propagation
        function propa (obj, state, eq_motion)
            global space time
            
            if(isempty(obj.frc_old)) % Initialization by Strang splitting
                
                obj.frc_old = state.frc;
                propa(tmp.traj.propagators.leap_frog, state, eq_motion);
                
            else
                
                for d = 1:space.n_dim
                    state.pos{d} = state.pos{d} + time.steps.s_delta * state.mom{d} / space.dof{d}.mass ...
                        + 1/6 * time.steps.s_delta^2 * (4 * state.frc{d} - obj.frc_old{d}) / space.dof{d}.mass;
                end
                
                frc_current = state.frc;
                eval_V_F (eq_motion, state, time.steps.s_delta, 0)
                
                for d = 1:space.n_dim
                    state.mom{d} = state.mom{d} + 1/6 * time.steps.s_delta * ...
                        (2 * state.frc{d} + 5 * frc_current{d} - obj.frc_old{d});
                end
                
                obj.frc_old = frc_current;
                
            end
            
        end
        
    end
end


