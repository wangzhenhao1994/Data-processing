%--------------------------------------------------------------------------
%
% Normalize (coupled) wavefunction(s)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-.... Burkhard Schmidt
%
% see the README file for license details.

function psi_norm = normalize ( psi_in )
global space hamilt

% Approximate normalization integral
norm2 = 0;
for m = 1:hamilt.coupling.n_eqs
    norm2 = norm2 + sum(abs(psi_in{m}(:)).^2 .* space.weight(:));
end

% Preallocate
psi_norm = cell(hamilt.coupling.n_eqs,1);

% Normalize each of the (coupled) wavefunction(s)
for m = 1:hamilt.coupling.n_eqs
    psi_norm{m} = psi_in{m} / sqrt(norm2);
end

end

