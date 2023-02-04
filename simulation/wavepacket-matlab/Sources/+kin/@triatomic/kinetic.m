% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function kinetic (obj, psi, new)

% We do three things:
%
% 1. Convert the wavefunction to FBR
% 2. Multiply with our kinetic energy grid
% 3. Convert it back to DVR
%
% If new is set to true, we apply everything to psi.new, otherwise
% to psi.dvr.
%
% Note that we transparently call kin_init if the grid is not yet set,
% so you might have to catch the initialised object in the output.

global space hamilt

% If kin_apply is called without previous initialisation, initialise it here.
if isempty (obj.grid)
	init_kin (obj, 1);
end


%% 1. Convert Wavefunction to FBR
for m = 1:hamilt.coupling.n_eqs
    if new
        psi.new{m} = dvr2fbr(space.dof{obj.dof(1)}, psi.new{m});
        psi.new{m} = dvr2fbr(space.dof{obj.dof(2)}, psi.new{m});
    else
        psi.new{m} = dvr2fbr(space.dof{obj.dof(1)}, psi.dvr{m});
        psi.new{m} = dvr2fbr(space.dof{obj.dof(2)}, psi.new{m});
    end
end

%% 2. Apply the grid
for m = 1:hamilt.coupling.n_eqs
    psi.new{m} = psi.new{m} .* obj.grid;
end

%% 3. Convert back
for m = 1:hamilt.coupling.n_eqs
    psi.new{m} = fbr2dvr(space.dof{obj.dof(1)}, psi.new{m});
    psi.new{m} = fbr2dvr(space.dof{obj.dof(2)}, psi.new{m});
end
