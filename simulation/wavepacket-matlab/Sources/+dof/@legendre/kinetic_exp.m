%--------------------------------------------------------------------------
% Performs a propagation step in the split operator method.
% Basically the same as grid_kinetic, but uses
% obj.kin_expo for the multiplication.
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function kinetic_exp( obj, psi )

global hamilt

if obj.nokin
    return
end

for m = 1:hamilt.coupling.n_eqs
    % 1. Expand in the FBR basis
    [psi.dvr{m}, permutation, shapedims] = obj.shape(psi.dvr{m});
    psi.dvr{m} = obj.trafo2fbr * psi.dvr{m};

    % 2. Apply the prepared propagator
    psi.dvr{m} = obj.kin_expo .* psi.dvr{m};

    % 3. Reconstruct the wave function in the pseudospectral basis.
    psi.dvr{m} = obj.trafo2dvr * psi.dvr{m};
    psi.dvr{m} = obj.shape_back(psi.dvr{m}, permutation, shapedims);
end
