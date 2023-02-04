%------------------------------------------------------------------------
% Obtain density matrix from 1-dimensional wavepacket propagation
%-------------------------------------------------------------------------
function wave_redu_1d (~,psi,m)
global space
    
for ii=1:space.dof{1}.n_pts
    for jj=1:space.dof{1}.n_pts
        psi.redu{m,1}(ii,jj) = conj(psi.dvr{m}(ii)) .* psi.dvr{m}(jj);
    end
end
