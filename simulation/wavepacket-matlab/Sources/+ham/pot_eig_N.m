%--------------------------------------------------------------------------
%
% Solve eigenproblem for the instantaneous (effective) potential matrix
%
% For the case of more than two coupled Schrödinger equations,
% the eigenvalue problem has be solved numerically
%
% Input argument "e" is the electric field (vector)
% Input argument "save_trafo" toggles saving of transformation matrices
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-20.. Burkhard Schmidt's group
%               2008,2010 Ulf Lorenz
%
% see the README file for license details.

function pot_eig_N ( e, save_trafo )
global hamilt space

% Pre-allocate storage for eigenvalues
hamilt.eig_val = cell (hamilt.coupling.n_eqs,1); % cell vector
for m = 1:hamilt.coupling.n_eqs
    hamilt.eig_val{m} = zeros(size(space.dvr{1}));
end

% Pre-allocate storage for eigenvectors
if save_trafo
    hamilt.eig_vec = cell (hamilt.coupling.n_eqs);   % cell matrix
    for m = 1:hamilt.coupling.n_eqs
        for n = 1:hamilt.coupling.n_eqs
            hamilt.eig_vec{m,n} = zeros(size(space.dvr{1}));
        end
    end
end

for gr=1:space.n_tot
    pot_mat = zeros(hamilt.coupling.n_eqs);
    
    % Extract diabatic potentials from cell arrays
    for m=1:hamilt.coupling.n_eqs
        if ~isempty ( hamilt.pot{m,m}.dvr ) 
            pot_mat(m,m) = hamilt.pot{m,m}.dvr(gr);
        end
            
        % Permanent dipole moments along p = x|y
        if isfield(hamilt, 'dip')
            for p = 1:length(hamilt.dip)
                if abs(e(p))>0
                    if ~isempty (hamilt.dip{p})
                        if ~isempty ( hamilt.dip{p}{m,m}.dvr )
                            pot_mat(m,m) = pot_mat(m,m) - e(p) * hamilt.dip{p}{m,m}.dvr(gr);
                        end
                    end
                end
            end
        end
        
        % Polarizabilities along p = x|y
        if isfield(hamilt, 'pol')
            for p = 1:size(hamilt.pol,1)
                if abs(e(p))>0
                    for q = p:size(hamilt.pol,2)
                        if abs(e(q))>0
                            if ~isempty (hamilt.pol{p,q}.dvr)
                                if ~isempty ( hamilt.pol{p,q}{m,m}.dvr )
                                    pot_mat(m,m) = pot_mat(m,m) - e(p)*e(q)/2 * hamilt.pol{p,q}{m,m}.dvr(gr);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % Three possible coupling mechanisms
        for n=m+1:hamilt.coupling.n_eqs
            
            % Diabatic potential coupling
            if ~isempty ( hamilt.pot{m,n}.dvr ) 
                pot_mat(m,n) = pot_mat(m,n) + hamilt.pot{m,n}.dvr(gr);
            end
            
            % Transition dipole moments along p = x|y
            if isfield(hamilt, 'dip')
                for p = 1:length(hamilt.dip)
                    if abs(e(p))>0
                        if ~isempty (hamilt.dip{p})
                            if  ~isempty(hamilt.dip{p}{m,n}.dvr)
                                pot_mat(m,n) = pot_mat(m,n) - e(p) * hamilt.dip{p}{m,n}.dvr(gr);
                            end
                        end
                    end
                end
            end
            
            % Polarizabilities along p = x|y
            if isfield(hamilt, 'pol')
                for p = 1:size(hamilt.pol,1)
                    if abs(e(p))>0
                        for q = p:size(hamilt.pol,2)
                            if abs(e(q))>0
                                if ~isempty (hamilt.pol{p,q})
                                    if  ~isempty(hamilt.pol{p,q}{m,n}.dvr)
                                        pot_mat(m,n) = pot_mat(m,n) - e(p)*e(q)/2 * hamilt.pol{p,q}{m,n}.dvr(gr);
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % Make potential matrix Hermitian
            pot_mat(n,m) = conj(pot_mat(m,n));
        end
    end
    
    % Diagonalize potential matrix obtained above
    % Columns of eig_vec are the right eigenvectors corresponding to eig_val
    [eig_vec, eig_val] = eig(pot_mat);
            
    % Save adiabatic potentials in cell vectors
    for m=1:hamilt.coupling.n_eqs
        hamilt.eig_val{m}(gr) = eig_val(m,m);
    end
    
    % Save adiabatic<=>diabatic transformation in cell arrays
    if save_trafo
        for m = 1:hamilt.coupling.n_eqs
            for n = 1:hamilt.coupling.n_eqs
                hamilt.eig_vec{m,n}(gr) = eig_vec (m,n);
            end
        end
    end
    
end
    
% Compensate sudden changes of eigenvector signs
if ~save_trafo
    return
end
old = zeros(hamilt.coupling.n_eqs,1);
new = zeros(hamilt.coupling.n_eqs,1);
arraySize = size(space.dvr{1});
switch space.n_dim
    
    case 1
        for i1 = 2:space.dof{1}.n_pts
            i_new = i1;
            i_old = i1-1;
            smoothen (i_new, i_old)
        end
        
    case 2
        
        % Along 1st coordinate
        for i2 = 1:space.dof{2}.n_pts
        for i1 = 2:space.dof{1}.n_pts
            i_new = sub2ind( arraySize, i1,   i2 );
            i_old = sub2ind( arraySize, i1-1, i2 );
            smoothen (i_new, i_old)
        end
        end
        
        % Along 2nd coordinate
        for i1 = 1:space.dof{1}.n_pts
        for i2 = 2:space.dof{2}.n_pts
            i_new = sub2ind( arraySize, i1, i2   );
            i_old = sub2ind( arraySize, i1, i2-1 );
            smoothen (i_new, i_old)
        end
        end
        
    case 3
        
        % Along 1st coordinate
        for i2 = 1:space.dof{2}.n_pts
        for i3 = 1:space.dof{3}.n_pts
        for i1 = 2:space.dof{1}.n_pts
            i_new = sub2ind( arraySize, i1,   i2, i3 );
            i_old = sub2ind( arraySize, i1-1, i2, i3 );
            smoothen (i_new, i_old)
        end
        end
        end
        
        % Along 2nd coordinate
        for i1 = 1:space.dof{1}.n_pts
        for i3 = 1:space.dof{3}.n_pts
        for i2 = 2:space.dof{2}.n_pts
            i_new = sub2ind( arraySize, i1, i2,   i3 );
            i_old = sub2ind( arraySize, i1, i2-1, i3 );
            smoothen (i_new, i_old)
        end
        end
        end
        
        % Along 3rd coordinate
        for i1 = 1:space.dof{1}.n_pts
        for i2 = 1:space.dof{2}.n_pts
        for i3 = 2:space.dof{3}.n_pts
            i_new = sub2ind( arraySize, i1, i2, i3   );
            i_old = sub2ind( arraySize, i1, i2, i3-1 );
            smoothen (i_new, i_old)
        end
        end
        end
        
    case 4
        
        % Along 1st coordinate
        for i2 = 1:space.dof{2}.n_pts
        for i3 = 1:space.dof{3}.n_pts
        for i4 = 1:space.dof{4}.n_pts
        for i1 = 2:space.dof{1}.n_pts
            i_new = sub2ind( arraySize, i1,   i2, i3, i4 );
            i_old = sub2ind( arraySize, i1-1, i2, i3, i4 );
            smoothen (i_new, i_old)
        end
        end
        end
        end
        
        % Along 2nd coordinate
        for i1 = 1:space.dof{1}.n_pts
        for i3 = 1:space.dof{3}.n_pts
        for i4 = 1:space.dof{4}.n_pts
        for i2 = 2:space.dof{2}.n_pts
            i_new = sub2ind( arraySize, i1, i2,   i3, i4 );
            i_old = sub2ind( arraySize, i1, i2-1, i3, i4 );
            smoothen (i_new, i_old)
        end
        end
        end
        end
        
        % Along 3rd coordinate
        for i1 = 1:space.dof{1}.n_pts
        for i2 = 1:space.dof{2}.n_pts
        for i4 = 1:space.dof{4}.n_pts
        for i3 = 2:space.dof{3}.n_pts
            i_new = sub2ind( arraySize, i1, i2, i3,   i4 );
            i_old = sub2ind( arraySize, i1, i2, i3-1, i4 );
            smoothen (i_new, i_old)
        end
        end
        end
        end
        
        % Along 4th coordinate
        for i1 = 1:space.dof{1}.n_pts
        for i2 = 1:space.dof{2}.n_pts
        for i3 = 1:space.dof{3}.n_pts
        for i4 = 2:space.dof{4}.n_pts
            i_new = sub2ind( arraySize, i1, i2, i3, i4   );
            i_old = sub2ind( arraySize, i1, i2, i3, i4-1 );
            smoothen (i_new, i_old)
        end
        end
        end
        end
        
    otherwise
        prt.error ('Smoothing of Dia2Adi transformation only up to 4 dimensions')
end


    function smoothen (ii_new, ii_old)
        for nn = 1:hamilt.coupling.n_eqs % loop over columns=eigenvectors
            for mm=1:hamilt.coupling.n_eqs % extract n-th column
                old(mm) = hamilt.eig_vec{mm,nn}(ii_old);
                new(mm) = hamilt.eig_vec{mm,nn}(ii_new);
            end
            if  norm ( old + new ) < norm ( old - new )
                for mm=1:hamilt.coupling.n_eqs % change sign of n-th column
                    hamilt.eig_vec{mm,nn}(ii_new) = - hamilt.eig_vec{mm,nn}(ii_new);
                end
            end
        end
        
    end


end

