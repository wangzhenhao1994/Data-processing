% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%               2012 Ulf Lorenz
%
% see the README file for license details.

function kinetic (obj, psi, new)

% We do three things:
%
% 1. Convert the wavefunction to a mixed FBR representation for the
%    Legendre grid and a DVR one for the other grids.
% 2. Multiply with our kinetic energy grid. Since it is diagonal, this
%    is done pointwise instead of a real matrix multiplication.
% 3. Convert it back to pure DVR
%
% If new is set to true, we apply everything to psi.new, otherwise
% to psi.dvr. The result is stored in psi.new.
%
% Note that we transparently call kin_init if the grid is not yet set,
% so you might have to catch the initialised object in the output.

global space hamilt

% If kin_apply is called without previous initialisation, initialise it here.
if isempty (obj.grid)
    init_kin (obj, 1);
end

for m = 1:hamilt.coupling.n_eqs
    % DVR => FBR
    if new
        fbr = dvr2fbr(space.dof{obj.dof_c}, psi.new{m});
    else
        fbr = dvr2fbr(space.dof{obj.dof_c}, psi.dvr{m});
    end

    % Apply the grid
    fbr = fbr .* obj.grid;

    % Convert it back.
    psi.new{m} = fbr2dvr( space.dof{obj.dof_c}, fbr );
end
