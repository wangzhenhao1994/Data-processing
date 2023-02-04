% The kinetic Hamiltonian is applied here in three steps
% 1. Transform the wave function to FBR
% 2. Multiply with a matrix representing the kinetic operator in FBR
% 3. Transform back the wave function to DVR
%
% Note that we use the internal shape and shape_back functions, this
% saves a few reshapes. The disadvantage is that obj.kin looks quite strange.

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function kinetic ( obj, psi, new )

global hamilt

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
    % 1. Transform the wave function to FBR
    if new
        [psi.new{m}, permutation, shapedims] = obj.shape(psi.new{m});
    else
        [psi.new{m}, permutation, shapedims] = obj.shape(psi.dvr{m});
    end
    psi.new{m} = obj.trafo2fbr * psi.new{m};

    % 2. Multiply with a matrix representing the kinetic operator in FBR
    psi.new{m} = obj.kin .* psi.new{m};

    % 3. Transform back the wave function to DVR
    psi.new{m} = obj.trafo2dvr * psi.new{m};
    psi.new{m} = obj.shape_back(psi.new{m}, permutation, shapedims);
end
