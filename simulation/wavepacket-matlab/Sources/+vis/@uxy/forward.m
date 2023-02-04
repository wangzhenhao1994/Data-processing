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

function forward( obj, first, last )
global control

figure(7)

% (1) Plot input (=control fields) u(t)
if isfield(control,'u')
    subplot (obj.numPlots,1,1);
    
    for len=1:control.u.dim
        h = plot(control.t.steps(first:last),control.u.forward(len,first:last));
        set(h, ...
            'LineStyle', obj.patterns{1}, ...
            'Color',     obj.colors(len,:), ...
            'LineWidth', obj.l_thin)
    end
    
    title (control.title);
    
end

% (2) Plot state vector x(t)
% Note: x(t) is shifted with respect to equilibrium x_e whereas z(t) is not
subplot (obj.numPlots,1,obj.numPlots-1);

for len=1:min(control.x.dim,control.x.max)
    h = plot (control.t.steps(first:last), abs(control.x.forward(len,first:last)+control.x.equilib(len)));
    set(h, ...
        'LineStyle', obj.patterns{1}, ...
        'Color',     obj.colors(len,:), ...
        'LineWidth', obj.l_thin)
end

if ~isfield(control,'u')
    title (control.title);
end

% (3) Plot output (=observables) y(t)
subplot (obj.numPlots,1,obj.numPlots);

for len=1:control.y.dim
    h = plot (control.t.steps(first:last), abs(control.y.forward(len,first:last)));
    set(h, ...
        'LineStyle', obj.patterns{1}, ...
        'Color',     obj.colors(len,:), ...
        'LineWidth', obj.l_thin)
end

% Update figure window
drawnow;

% Optionally add a frame to the movie file
if control.plot.mov
    frame = getframe(gcf);
    writeVideo(obj.writerObj,frame);
end



end

