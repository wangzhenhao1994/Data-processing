%-------------------------------------------------------------------------
%
% Propagate objects of class traj (trajectory bundles) 
% by one substep of size time.steps.s_delta
%
% Yoshida (fourth order):
% -----------------------
%
%     Reference: H. Yoshida
%     Phys. Lett. A 150, 262 (1990)
%     DOI:10.1016/0375-9601(90)90092-3
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

classdef yoshida < tmp.traj.propagators.generic & handle
    
    properties (Access = private)
        c           % coefficients
        d           % coefficients
    end 
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = yoshida
            obj.string = 'Propagator: Yoshida''s algorithm (fourth order)';
        end
        
        % Initialize QC propagator
        function init (obj, state, eq_motion)
            
            % Inherit from "generic" super class
            init@tmp.traj.propagators.generic(obj, state, eq_motion);
            
            obj.c = zeros(4,1);
            obj.d = zeros(4,1);
            
            obj.c([1,4]) = 1 / (2 * ( 2 - 2^(1/3) ) );
            obj.c([2,3]) = ( 1 - 2^(1/3) ) / (2 * ( 2 - 2^(1/3) ) );
            
            obj.d([1,3])  = 1 / ( 2 - 2^(1/3) ) ;
            obj.d(2)      = - 2^(1/3) / ( 2 - 2^(1/3) ) ;
            obj.d(4)      = 0;      
            
        end

        % Perform QC propagation
        function propa (obj, state, eq_motion)
            global space time
            
            for indices = 1:4
                
                % Propagate positions
                for d_d = 1:space.n_dim
                    state.pos{d_d} = state.pos{d_d} + obj.d(indices) .* time.steps.s_delta * state.mom{d_d} / space.dof{d_d}.mass;
                end
                
                % Get potentials and forces
                eval_V_F (eq_motion, state, obj.d(indices) .* time.steps.s_delta, 0)
                
                % Propagate momenta
                for d_d = 1:space.n_dim
                    state.mom{d_d} = state.mom{d_d} + obj.c(indices) .* time.steps.s_delta * state.frc{d_d};
                end
            end
            
        end
        
    end
end