%-------------------------------------------------------------------------
%
% Solve the time-independent Schroedinger equation to get eigenstates
% and energies in position representation by using DVR/FBR techniques 
%
% Part 1/3: Set up the Hamiltonian matrix in DVR space
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function setup (obj)

global hamilt space

prt.disp (' ')
prt.disp('Start setting up matrix ...')
tim = cputime;

% Initialize matrix: Full or sparse
mat_size = space.n_tot * hamilt.coupling.n_eqs;
if obj.storage == 'f'
    obj.matrix = zeros(mat_size);
else
    obj.matrix = sparse(mat_size, mat_size);
end

% Total number (N) of grid points should be even!
N = space.n_tot;

% Create kinetic energy matrix from both the grid's internal kinetic energy
% operators and external kinetic energy operators. If we solve the TISE for
% coupled equations, the resulting matrix should have a blockwise-diagonal form;
% we fill only the left upper square and copy the result to the other blocks.

for k = 1:space.n_dim
    prt.disp(['kinetic energy for dof : ' num2str(k)])
    obj.matrix(1:N, 1:N) = ...
        obj.matrix(1:N,1:N) + kinetic2dvr(space.dof{k});
end

if isfield(hamilt, 'kin')
    for ii = 1:length(hamilt.kin)
        prt.disp(['custom kinetic energy operator : ' num2str(ii)]);
        obj.matrix(1:N, 1:N) = ...
            obj.matrix(1:N,1:N) + kinetic2dvr(hamilt.kin{ii});
    end
end

for m = 2:hamilt.coupling.n_eqs
    obj.matrix((m-1)*N+1:m*N, (m-1)*N+1:m*N) = ...
        obj.matrix(1:N, 1:N);
end

% Add the potential energy and the diabatic coupling elements
% Since the potential energy is not given as a diagonal matrix, but as a grid, we have
% to make a matrix out of it.
prt.disp('potential energy')

for m = 1:hamilt.coupling.n_eqs
    for n = m:hamilt.coupling.n_eqs
        if isempty(hamilt.pot{m,n}.dvr)
            continue;
        end

        pot = reshape(hamilt.pot{m,n}.dvr, N, 1); % or hamilt.pot{m,n}.dvr(:) !?!
        if obj.storage == 's'
            pot = sparse(pot);
        end

        obj.matrix((m-1)*N+1:m*N, (n-1)*N+1:n*N) = ...
            obj.matrix((m-1)*N+1:m*N, (n-1)*N+1:n*N) ...
            + diag(pot);

        if m ~= n
            % diabatic couplings should be symmetric
            obj.matrix((n-1)*N+1:n*N, (m-1)*N+1:m*N) = ...
                diag(pot);
        end
    end
end

% Cut-off
obj.matrix(abs(obj.matrix) < obj.cutoff) = 0;

prt.disp (['Finished after [CPU seconds] : ' num2str(cputime-tim)])
prt.disp (' ')

