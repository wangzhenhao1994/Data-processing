%----------------------------------------------------------
% Plot color mapping
%----------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function show_colors (obj)

x = linspace ( -1, 1, 50 );
y = ones (size(x));
c = linspace ( 0, 1, 50 );
vis.styles.plot_color ( x, y, c, [1 0 0], obj.l_extra, 0, 0)
axis ( [-1 1 0.9 1.1] )
set(gca, ...
    'XTick',          -0.5:0.5:1, ...
    'YTick',         [-123456789 123456789], ... % Fake: suppress tick labels
    'YAxisLocation', 'right', ...
    'LineWidth',     obj.l_thick, ...
    'FontName',      obj.f_name, ...
    'FontSize',      obj.f_large, ...
    'FontWeight',    obj.f_heavy)
xlabel ('phase [\pi]')
ylabel ('color map')

end
