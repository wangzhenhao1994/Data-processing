%--------------------------------------------------------------------------
%
% Plot evolution of input/state/output of bilinear control problem 
% in one figure with 2 or 3 subplots 
% (1) input=control field(s) u(t) - not drawn for the field-free case - 
% (2) state vector x(t)
% (3) output=observable(s) y(t) 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function equilib ( obj )

global control time

figure(7);

% Clear figure and draw logos in all four corners
clf;
show_logo (obj);

% (1) Input (=control fields) u(t)
if isfield(control,'u')
    subplot(obj.numPlots,1,1);
    
    axis ( [[control.t.steps(1) control.t.steps(end)] [-1 1]*2*abs(time.efield.max_ampli) ] )
    hold on
    for len=1:control.u.dim
        h = line([control.t.steps(1) control.t.steps(end)],[0 0]);
        set(h, ...
            'LineStyle',   ':', ...
            'Color',       obj.colors(len,:), ...
            'LineWidth',   obj.l_thin, ...
            'DisplayName', control.u.legends{len})
    end
    
    set(gca, ...
        'LineWidth',  obj.l_thick, ...
        'FontName',   obj.f_name, ...
        'FontSize',   obj.f_large, ...
        'FontWeight', obj.f_heavy)

    xlabel('time')
    ylabel('control field(s) u(t)')
    box on
    if obj.legends
        legend('Location','NorthWest')
    end
    
end

% (2) State vector x(t)
subplot(obj.numPlots,1,obj.numPlots-1);

axis ( [ [control.t.steps(1) control.t.steps(end)] [-0.1 1.1]] )
hold on
for len=1:min(control.x.dim,control.x.max)
    h = line([control.t.steps(1) control.t.steps(end)],real(control.x.equilib(len))*[1 1]);
    set(h, ...
        'LineStyle', '--', ...
        'Color',     obj.colors(len,:), ...
        'LineWidth', obj.l_thin, ...
        'DisplayName', control.x.legends{len})
end

set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)

xlabel('time')
ylabel('state vector |x(t)|')
box on
if obj.legends
    legend('Location','NorthWest')
end

% (3) Output (=observables) y(t)
subplot (obj.numPlots,1,obj.numPlots);

axis ( [ [control.t.steps(1) control.t.steps(end)] [-0.1 1.1]] )
hold on
for len=1:control.y.dim
    h = line([control.t.steps(1) control.t.steps(end)],real(control.y.equilib(len))*[1 1]);
    set(h, ...
        'LineStyle',   '--', ...
        'Color',       obj.colors(len,:), ...
        'LineWidth',   obj.l_thin, ...
        'DisplayName', control.y.legends{len})
end

set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)

xlabel('time')
ylabel('observable(s) |y(t)|')
box on
if obj.legends
    legend('Location','NorthWest')
end


end

