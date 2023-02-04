%--------------------------------------------------------------------
% Plot external electric field vs. time
%--------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function efield (obj)
global time

if time.efield.dressed
    factor = 2;
else
    factor = 1;
end

% Set range of time (sub!)steps
if time.steps.offset == 1
    steps = 1; % first main time step
else
    steps = 1:time.steps.offset+time.steps.s_number; % all later main time steps
end

% Plot x-component, y-component of electric field
if time.efield.has_x
    h = plot ( time.steps.s_grid(steps), factor * real(time.efield.grid{1}(steps)) );
    set(h, ...
        'LineStyle', obj.patterns{2}, ...
        'Color',     'black', ...
        'LineWidth', obj.l_thick)
    hold on
end
if time.efield.has_y
    h = plot ( time.steps.s_grid(steps), factor * real(time.efield.grid{2}(steps)) );
    set(h, ...
        'LineStyle', obj.patterns{3}, ...
        'Color',     'black', ...
        'LineWidth', obj.l_thick)
    hold off
end

% Legend explaining the line styles
if obj.legends && time.efield.has_x && time.efield.has_y
    legend('F_x','F_y','Location','SouthWest')
end

% Axes and labels
if time.efield.dressed
    axis ( [ 0 time.steps.t_total  [-1.1 +1.1]*time.efield.max_ampli ] )
else
    line ([0 time.steps.t_total],[0 0],'Color','black','LineWidth', obj.l_thin)
    axis ( [ 0 time.steps.t_total -abs(time.efield.max_ampli)*1.1 +abs(time.efield.max_ampli)*1.1 ] )
end

set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ... 
    'FontWeight', obj.f_heavy)
xlabel ('t')
ylabel ('F(t)' )

end
        
