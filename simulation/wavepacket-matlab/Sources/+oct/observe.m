%--------------------------------------------------------------------------
%
% Calculates values of observables
% which can be linear (state.C) and/or quadratic (state.D) 
%
% y = C * x + x * D * x
%
% Depending on variable 'action' the following is achieved:
% 'initial'  : initialization 
% 'forward'  : forward propagation 
% 'backward' : backward propagation 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2012 Burkhard Schmidt, Jeremy Rodriguez, Ulf Lorenz
%
% see the README file for license details.

function observe (action, step)

global control state;

switch lower(action)
    
    %% Only upon first call: get/check number of observables
    case 'initial'
        
        % At least one of C or D should exist
        if isempty(state.C) && isempty(state.D)
            prt.error ('either cell vector C or D need to exist')
        end
        
        % However, C and D shouldn't both exist (to be generalized later ...)
        if ~isempty(state.C) && ~isempty(state.D)
            prt.error ('cell vectors C and D should not both exist')
        end
        
        % Lengths of cell vectors C or D (the latter should be Hermitian)
        if ~isempty(state.C)
            state.len_CD = length(state.C);
        end
        if ~isempty(state.D)
            state.len_CD = length(state.D);
            for len=1:state.len_CD
%                 if ~ishermitian (state.D{len})
%                     prt.error ('All observable matrices D should be Hermitian')
%                 end
            end
        end
        
        %% Calculate observables from C and/or D from state vector x(t)
    case 'forward'
        % Note: x(t) is shifted with respect to equilibrium x_e whereas z(t) is not
        x = control.x.forward (:,step)+control.x.equilib;        
        y = zeros (state.len_CD,1);
        for len = 1:state.len_CD
            
            if ~isempty(state.C) % linear observable: <c|x> may be complex
                if ~state.Q{len}
                    y(len) = y(len) + real(state.C{len}*x);
                else
                    y(len) = y(len) +  abs(state.C{len}*x)^2;
                end
            end
            if ~isempty(state.D) % quadratic observable: real if D is Hermitian
                y(len) = y(len) + dot( x, state.D{len} * x);
            end
        end
        control.y.forward(:,step) = y;
        
        %% Calculate observables from C and/or D from Lagrange multiplier z(t)
    case 'backward'
        % Note: x(t) is shifted with respect to equilibrium x_e whereas z(t) is not
        z = control.x.backward(:,step);
        y = zeros (state.len_CD,1);
        for len = 1:state.len_CD
            if ~isempty(state.C) % linear observable: <c|z> may be complex
                if ~state.Q{len}
                    y(len) = y(len) + real(state.C{len}*z);
                else
                    y(len) = y(len) +  abs(state.C{len}*z)^2;
                end
            end
            if ~isempty(state.D) % quadratic observable: real if D is Hermitian
                y(len) = y(len) + dot( z, state.D{len} * z);
            end
        end
        control.y.backward(:,step) = y;
        
    otherwise
        prt.error('Calling observe.m with invalid keyword')
end

end
