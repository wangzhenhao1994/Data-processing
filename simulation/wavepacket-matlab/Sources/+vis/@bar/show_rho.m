%-------------------------------------------------
% Plot "rho" matrices as 3-D bar graphs
%-------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-.... Burkhard Schmidt
%
% see the README file for license details.

function show_rho ( obj, state, ~ )
global control info 

%% Prepare data
n = sqrt(length(state.x));
x = state.x+state.x_equilib;
if isfield (control.lvne,'order') && strcmpi(control.lvne.order,'df')
    U=math.cw2df(n);
    x = U' * x;
end
x = reshape(x,[n,n]);

switch lower(obj.complex)
    case 'abs2'
        b = abs (x);                 % not squaring here
        s = '|\rho_{nm}|';
    case 'real'
        b = real (x);
        s = 'Re(\rho_{nm})';
    case 'imag'
        b = imag (x);
        s = 'Im(\rho_{nm})';
    otherwise
        prt.error('Wrong choice for obj.complex')
end

%% Create 3-D bar graph and set view point in terms of azimuth and elevation
bar3 (b,'r')
view (obj.bar_view(1),obj.bar_view(2))


%% Axes, labels, etc
axis ([-1.0 n+1 -1.0 n+1 -0.1 +1.1])

ax = gca;
set (ax, 'LineWidth',  obj.l_thick);
set (ax, 'FontName',   obj.f_name);
set (ax, 'FontSize',   obj.f_large);
set (ax, 'FontWeight', obj.f_heavy);

title ( {info.header1;info.header2} )
xlabel ( 'n' )
ylabel ( 'm' )
zlabel ( s )

end
