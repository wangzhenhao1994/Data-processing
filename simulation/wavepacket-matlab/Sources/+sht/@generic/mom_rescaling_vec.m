%--------------------------------------------------------------------------
%
% Rescaling momenta of classical particles to enforce
% energy conservation in surface hopping trajectories
%
% For a transition from state m to state n, try
% momentum rescaling along a cell-vector "vec" 
% Note that mom, D and vec should all be based
% on the same number of trajectories. 
%
% Note that upon hopping from lower to upper states
% it may be that the momenta are not sufficient in
% which case the energy cannot be conserved, hence
% the hopping is not performed ("frustrated hops")
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function [mom_new, allowed, forbidden] = mom_rescaling_vec (~,mom,D,vec,m,n)
global space

% Pre-allocate memory
n_hop   = length(mom{1});
mom_new = cell(size(mom));

% Potential energy before|after hopping
pot_m = D(m,:)';
pot_n = D(n,:)';

A = zeros (n_hop, 1);
B = zeros (n_hop, 1);

for d = 1:space.n_dim
    A = A + 0.5 / space.dof{d}.mass * vec{d}.^2;
    B = B + 1.0 / space.dof{d}.mass * vec{d}.* mom{d};
end

C = pot_n - pot_m;
D = B ./ (2*A);
E = C ./ A;

allowed   = (D.^2 - E >= 0);
forbidden = (D.^2 - E  < 0);

if ~isempty(allowed)
    ratio = - D + (-1).^(D<0) .* sqrt( (D.^2 - E) );
    for d = 1:space.n_dim
        mom_new{d}(allowed) = mom{d}(allowed) + ratio(allowed) .* vec{d}(allowed);
    end
end
end
