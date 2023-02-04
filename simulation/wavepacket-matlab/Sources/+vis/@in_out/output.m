%------------------------------------------------------------------------------------
% Output of a bilinear control system: observables versus time
%------------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-20.. Burkhard Schmidt
%
% see the README file for license details.

function output (obj,step)
global info state time

for len=1:size(state.y,2)
    h = plot (time.steps.m_grid(1:step), state.y(1:step,len));
    set(h, ...
        'LineStyle', obj.patterns{1}, ...
        'Color',     obj.colors(len,:), ...
        'LineWidth', obj.l_thick, ...
        'DisplayName', state.y_label{len})
    if len==1
        hold on
    end
end
hold off

% Set axis limits 
if ~isempty(obj.y_min)
    y_min = obj.y_min;
else
    y_min = -inf;
end

if ~isempty(obj.y_max)
    y_max = obj.y_max;
else
    y_max = +inf;
end

axis ( [ 0 time.steps.t_total y_min y_max ] )

% Legend explaining the curve colors
if obj.legends
    legend('Location','SouthEast') 
end

% Fonts and line styles
set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)

% Axis labels
switch lower (info.program)
    case {'qm_propa'} % Time dependent simulations
        xlabel ('t')
    case{'qm_bound'} % Time independent: Loop over eigenstates
        xlabel ('n')
    otherwise
        prt.error ('Wrong choice of program')
end
ylabel ( 'output: observables' )

% Plot title
if ~isempty(state.title)
    title([state.title int2str(length(state.x)) ' coupled ODEs'])
end

end
