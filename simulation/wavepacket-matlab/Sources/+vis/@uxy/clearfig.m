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

function clearfig( obj )
global control
        
figure(7)

% input (=control fields) u(t)
if isfield(control,'u')
    subplot (obj.numPlots,1,1);
    hold off
end

% state vector x(t)
subplot (obj.numPlots,1,obj.numPlots-1);
hold off

% output (=observables) y(t)
subplot (obj.numPlots,1,obj.numPlots);
hold off

end

