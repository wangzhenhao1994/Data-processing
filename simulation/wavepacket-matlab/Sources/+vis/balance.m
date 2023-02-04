%--------------------------------------------------------------------------
% Plot details of balancing transformation: 
% (1) Hankel singular values (HSVs)
%     Sigma:  Hankel singular values (HSVs) from transformation of WC, WO
% (2) Transformation matrix
%--------------------------------------------------------------------------

function balance (Sigma, WC, WO, S, T)

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Boris Schaefer-Bung, Burkhard Schmidt, Ulf Lorenz
%
% see the README file for license details.

global reduce state

figure(9);
clf;
thisplot = vis.styles;

% First subplot: Hankel singular values
subplot (1,2,1);
show_logo (thisplot);
set ( gca, ...
    'LineWidth',     thisplot.l_thick, ...
    'FontName',      thisplot.f_name,  ...
    'FontSize',      thisplot.f_large, ...
    'FontWeight',    thisplot.f_heavy )
semilogy(real(Sigma),                               'xk','Markersize',5)
hold on
semilogy(sqrt(abs(eig(WC*WO))),                     'o','Markersize',5)
semilogy(abs(real(diag(S*WC*S'))),'x','Markersize',5)
semilogy(abs(real(diag(T'*WO*T))),'o','Markersize',5)
legend(...
    'Hankel singular values', ...
    'sqrt(|eig(W_C*W_O)|)', ...
    'S*W_C*S''', ...
    'T''*W_O*T', ...
    'Location', 'NorthEast')
set ( gca, ...
    'LineWidth',     thisplot.l_thick, ...
    'FontName',      thisplot.f_name,  ...
    'FontSize',      thisplot.f_large, ...
    'FontWeight',    thisplot.f_heavy )
xlabel 'HSV index'
ylabel 'singular values \Sigma_i'
title({['Hankel Singular Values \Sigma: ' upper(reduce.balance.transform) '-method'] state.title})
axis([0 size(S,1) 1e-20 1.01])
hold off

% Second subplot: balancing transformation
subplot (1,2,2);
image(log10(abs(S')),'CDataMapping','scaled')
set(gca,'YDir','normal')
colorbar ('eastoutside')
set ( gca, ...
    'LineWidth',     thisplot.l_thick, ...
    'FontName',      thisplot.f_name,  ...
    'FontSize',      thisplot.f_large, ...
    'FontWeight',    thisplot.f_heavy )
title ('Balancing transform: log_{10}(abs(S^T))')
xlabel 'singular vectors'
ylabel 'original basis'

end
