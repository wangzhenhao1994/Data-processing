%-------------------------------------------------
% Plot "ket" vectors as 2-D bar graphs
%-------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-.... Burkhard Schmidt
%
% see the README file for license details.

function show_ket ( obj, state, ~ )
global info 

%% Prepare data
n = length(state.x);
x = state.x+state.x_equilib;
switch lower(obj.complex)
    case 'abs2'
        b = abs (x) .^2;
        s = '|\psi_n|^2';
    case 'real'
        b = real (x);
        s = 'Re(\psi_n)';
    case 'imag'
        b = imag (x);
        s = 'Im(\psi_n)';
end

%% Create 2-D bar graphs
bar (0:n-1,b)

%% Axes, labels, etc
axis ([-1.0 length(state.x) -0.1 +1.1])

ax = gca;
set (ax, 'LineWidth',  obj.l_thick);
set (ax, 'FontName',   obj.f_name);
set (ax, 'FontSize',   obj.f_large);
set (ax, 'FontWeight', obj.f_heavy);

title ( {info.header1;info.header2} )
xlabel ( 'n' )
ylabel ( s )

end

