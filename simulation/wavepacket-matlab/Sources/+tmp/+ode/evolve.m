%--------------------------------------------------------------------------
% Right-hand-side of bilinear control equation for state vector x(t):
%
%  d                     m                         
%  -- x(t) = A x(t) + i Sum u (t) ( N x(t) + b )           
%  dt                   k=1   k       k        k
%
% where A is used to describe the dynamics of the unperturbed system and 
% where the external control field(s) u_k(t) interact through N_k and b_k
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Boris Schaefer-Bung, Burkhard Schmidt, Ulf Lorenz
%               2012 Burkhard Schmidt, Jeremy Rodriguez
%
% see the README file for license details.

function x_dot = evolve (t, x)

global state time

% dynamics of the unperturbed system
x_dot = state.A*x;

% Get external control field
if isfield (time,'pulse')
    u = efi.eval(t);
    
    % interaction with external control fields along x and or y
    for d=1:length(u)
        if abs(u{d})~=0
            x_dot = x_dot + 1i * u{d}  * ( state.N{d} * x  + state.B{d} );
        end
    end
    
end

