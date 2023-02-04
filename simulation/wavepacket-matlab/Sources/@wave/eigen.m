%--------------------------------------------------------------------------
%
% Extract eigenfunctions from selected columns(!) of eigenvector matrix 
% of the Hamiltonian and normalize
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%
% see the README file for license details.

function eigen ( psi, step )
global hamilt space

index = step + hamilt.eigen.start;

% With symmetry
if ~isempty(hamilt.eigen.symmetry) && ~strcmpi(hamilt.eigen.symmetry,'n')
    bloated = hamilt.eigen.transform' * hamilt.eigen.eig_vecs(:, index);
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr{m} = bloated((m-1)*space.n_tot+1:m*space.n_tot);
    end

% General case: without symmetry
else
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr{m} = hamilt.eigen.eig_vecs((m-1)*space.n_tot+1:m*space.n_tot,index);
    end
end

% Reshape eigenvector for >1 dimension
dims = zeros (1,space.n_dim); % row vector
if space.n_dim > 1
    for k = 1:space.n_dim
        dims(k) = space.dof{k}.n_pts;
    end
    
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr{m} = reshape(psi.dvr{m}, dims);
    end
end


% Normalize eigenfunction
psi.dvr = wave.normalize(psi.dvr);

end
