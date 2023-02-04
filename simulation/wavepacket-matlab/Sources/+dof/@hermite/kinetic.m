% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.


function kinetic ( obj, psi, new)

global hamilt

% The kinetic Hamiltonian is applied here in a single step, as the
% transformation to and from FBR has already been included in the
% definition of the internal matrix.

if isempty(obj.kin)
    init_kin(obj, 1);
end

if obj.nokin
    for m = 1:hamilt.coupling.n_eqs
        psi.new{m} = zeros(size(psi.dvr{m}));
    end
    return
end

for m = 1:hamilt.coupling.n_eqs
    % Shape the grid so that our DOF is first
    if new
        [psi.new{m}, permutation, shapedims] = obj.shape(psi.new{m});
    else
        [psi.new{m}, permutation, shapedims] = obj.shape(psi.dvr{m});
    end

    % Apply the kinetic energy operator
    psi.new{m} = obj.kin * psi.new{m};

    % Reshape the grid to the original form
    psi.new{m} = obj.shape_back(psi.new{m}, permutation, shapedims);
end
