%------------------------------------------------------------------------
% Reduce total density matrix (from 4-dimensional wavepacket propagation)
% to k-th dimension (k=1,2,3,4) by tracing over all the remaining dimensions
%-------------------------------------------------------------------------
function wave_redu_4d (~,psi,m)
global space

% First coordinate
for ii=1:space.dof{1}.n_pts
    for jj=1:space.dof{1}.n_pts
        psi.redu{m,1}(ii,jj) = sum ( sum ( sum ( conj(psi.dvr{m}(ii,:,:,:)) .* psi.dvr{m}(jj,:,:,:) ) ) );
    end
end

% Second coordinate
for ii=1:space.dof{2}.n_pts
    for jj=1:space.dof{2}.n_pts
        psi.redu{m,2}(ii,jj) = sum ( sum ( sum ( conj(psi.dvr{m}(:,ii,:,:)) .* psi.dvr{m}(:,jj,:,:) ) ) );
    end
end

% Third coordinate
for ii=1:space.dof{3}.n_pts
    for jj=1:space.dof{3}.n_pts
        psi.redu{m,3}(ii,jj) = sum ( sum ( sum ( conj(psi.dvr{m}(:,:,ii,:)) .* psi.dvr{m}(:,:,jj,:) ) ) );
    end
end

% Fourth coordinate
for ii=1:space.dof{4}.n_pts
    for jj=1:space.dof{4}.n_pts
        psi.redu{m,4}(ii,jj) = sum ( sum ( sum ( conj(psi.dvr{m}(:,:,:,ii)) .* psi.dvr{m}(:,:,:,jj) ) ) );
    end
end
