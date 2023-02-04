%-------------------------------------------------------------------------
%
% Transform (vector-valued) wavefunction from adiabatic to 
% diabatic representation. Within WavePacket, this routine 
% is used _only_ for the purpose of transforming an initial 
% wavefunction in an adiabatic representation to the diabatic
% representation for the first time. If you wish to transform
% the wave function later, use wave.adiabatic instead.
% 
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% see the README file for license details.

function diabatic ( psi )
global hamilt space

% Solve eigenproblem
switch hamilt.coupling.n_eqs
    case 1
        return
        
    case 2 % Analytic diagonalization for two coupled channels
        ham.pot_eig_2 ([0 0], 1)
        
    otherwise % Numerical diagonalization for more than two coupled channels
        ham.pot_eig_N ([0 0], 1)
        
end

% Save adiabatic wave functions 
for m = 1:hamilt.coupling.n_eqs
    psi.adi{m} = psi.dvr{m};
    psi.dvr{m} = zeros (size(space.dvr{1}));
end

% Transform to diabatic representation using eigenvector matrix
for m = 1:hamilt.coupling.n_eqs
    for n = 1:hamilt.coupling.n_eqs
        psi.dvr{m} = psi.dvr{m} + hamilt.eig_vec{m,n}.*psi.adi{n};
    end
end


end
