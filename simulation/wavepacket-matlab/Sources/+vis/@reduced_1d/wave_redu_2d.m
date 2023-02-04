%------------------------------------------------------------------------
% Reduce total density matrix (from 2-dimensional wavepacket propagation)
% to k-th dimension (k=1,2) by tracing over all the remaining dimensions
%-------------------------------------------------------------------------
function wave_redu_2d (~,psi,m)   
global space
    
% First coordinate
for ii=1:space.dof{1}.n_pts
    for jj=1:space.dof{1}.n_pts
        psi.redu{m,1}(ii,jj) = sum ( conj(psi.dvr{m}(ii,:)) .* psi.dvr{m}(jj,:) );
    end
end

% Second coordinate
for ii=1:space.dof{2}.n_pts
    for jj=1:space.dof{2}.n_pts
        psi.redu{m,2}(ii,jj) = sum ( conj(psi.dvr{m}(:,ii)) .* psi.dvr{m}(:,jj) );
    end
end
