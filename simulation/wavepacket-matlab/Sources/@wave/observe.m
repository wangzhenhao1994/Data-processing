%--------------------------------------------------------------------------
%
% For given wavefunctions (psi) and their representions (space) and for
% given a Hamiltonian operator (hamilt), this function calculates the
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

function observe ( obj, step )
global expect hamilt space

% Populations
wave_pop (expect.pop, obj, step);
wave_tot (expect.pop, obj, step);

if isfield(hamilt, 'amo')
    for k = 1:length(hamilt.amo)
        if ~isempty (hamilt.amo{k})
            wave_cha (expect.amo{k}, obj, step);
            wave_tot (expect.amo{k}, obj, step);
        end
    end
end

% Positions and momenta for each spatial dimension
for k = 1:space.n_dim
    wave_cha (expect.pos{k}, obj, step);
    wave_tot (expect.pos{k}, obj, step);
    wave_fbr (expect.mom{k}, obj, step);
    wave_tot (expect.mom{k}, obj, step);
end

% Potential and kinetic energy
wave_cha (expect.pot, obj, step);
wave_tot (expect.pot, obj, step);
wave_kin (expect.kin, obj, step);
wave_tot (expect.kin, obj, step);

% Total energy including couplings: Should be a real number
apply_ham (obj, [0 0],0); 
energy = 0;
for m = 1:hamilt.coupling.n_eqs
    energy = energy + sum ( conj(obj.dvr{m}(:)) .* obj.new{m}(:) .* space.weight(:) );
end
energy = energy / expect.pop.tot(step);
expect.total(step) = math.real(energy);
