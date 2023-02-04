%--------------------------------------------------------------------------
% Reduce total density (from 6-dimensional wavepacket propagation) to
% dimensions (1,2), (3,4), (5,6) by tracing over remaining dimensions
%--------------------------------------------------------------------------
function wave_redu_6d (~,state,m)
global space

% First and second coordinates
for ii=1:space.dof{1}.n_pts
    for jj=1:space.dof{2}.n_pts
        state.redu{m,1}(ii,jj) = sum ( sum ( sum ( sum ( abs(state.dvr{m}(ii,jj,:,:,:,:)))).^2  ) );
    end
end

% Third and fourth coordinates
for ii=1:space.dof{3}.n_pts
    for jj=1:space.dof{4}.n_pts
        state.redu{m,2}(ii,jj) = sum ( sum ( sum ( sum ( abs(state.dvr{m}(:,:,ii,jj,:,:)))).^2  ) );
    end
end

% Fifth and sixth coordinates
for ii=1:space.dof{5}.n_pts
    for jj=1:space.dof{6}.n_pts
        state.redu{m,3}(ii,jj) = sum ( sum ( sum ( sum ( abs(state.dvr{m}(:,:,:,:,ii,jj)))).^2  ) );
    end
end