%------------------------------------------------------------------------
% Obtain position densities from 2-dimensional wavepacket propagation
%-------------------------------------------------------------------------
function wave_redu_2d (~,state,m)   

state.redu{m,1} = abs(state.dvr{m}).^2;