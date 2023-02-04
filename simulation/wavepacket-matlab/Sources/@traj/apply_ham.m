%---------------------------------------------------------------------
%
% Application of Hamiltonian function to trajectory objects
%
%---------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20xy Burkhard Schmidt's group
%
% see the README file for license details.

function apply_ham ( obj )
global space

% Evaluate potential energy
% eval_V_F ( obj )

% Evaluate kinetic energy
obj.kin = zeros ( size ( obj.pos{1} ) );
for d = 1: space.n_dim
    eval_kin (  space.dof{d}, obj ) 
end

end

