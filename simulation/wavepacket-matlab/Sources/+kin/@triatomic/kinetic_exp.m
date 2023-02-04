% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function kinetic_exp ( obj, psi )

% We do three things:
%
% 1. Convert the wavefunction to FBR
% 2. Multiply with our exponentiated kinetic energy grid
% 3. Convert it back to DVR
%
% If new is set to true, we apply everything to psi.new, otherwise
% to psi.dvr.

global space hamilt


%% 1. Convert Wavefunction to FBR
for m = 1:hamilt.coupling.n_eqs
        psi.dvr{m} = dvr2fbr(space.dof{obj.dof(1)}, psi.dvr{m});
        psi.dvr{m} = dvr2fbr(space.dof{obj.dof(2)}, psi.dvr{m});
end

%% 2. Apply the grid
for m = 1:hamilt.coupling.n_eqs
    psi.dvr{m} = psi.dvr{m} .* obj.grid_exp;
end

%% 3. Convert back
for m = 1:hamilt.coupling.n_eqs
    psi.dvr{m} = fbr2dvr(space.dof{obj.dof(1)}, psi.dvr{m});
    psi.dvr{m} = fbr2dvr(space.dof{obj.dof(2)}, psi.dvr{m});
end
