%--------------------------------------------------------------------------
%
% For a given trajectory swarm (object of class "traj") and for a
% given Hamiltonian operator (hamilt), this function calculates the
% expectation values and uncertainties of all relevant observables
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function observe ( state, step )
global expect hamilt space

% Populations
for m=1:hamilt.coupling.n_eqs
    expect.pop.cha{m}(step) = nnz ( state.cha==m ) / state.n_p;
end
traj_tot (expect.pop, state, step);

if isfield(hamilt, 'amo')
    for k = 1:length(hamilt.amo)
        if ~isempty (hamilt.amo{k})
            traj_cha (expect.amo{k}, state, step);
            traj_tot (expect.amo{k}, state, step);
        end
    end
end

% Position and momenta for each spatial dimension
for k = 1:space.n_dim
     traj_cha (expect.pos{k}, state, step);
     traj_tot (expect.pos{k}, state, step);
     traj_cha (expect.mom{k}, state, step);
     traj_tot (expect.mom{k}, state, step);
end

% Potential and kinetic energy
traj_cha (expect.pot, state, step);
traj_tot (expect.pot, state, step);
traj_cha (expect.kin, state, step);
traj_tot (expect.kin, state, step);

% Total energy, not yet including couplings
energy = 0;
for m = 1:hamilt.coupling.n_eqs
    energy = energy + (expect.pot.cha{m}(step)+expect.kin.cha{m}(step)) * expect.pop.cha{m}(step);
end
energy = energy / expect.pop.tot(step);
expect.total(step) = energy;
