%--------------------------------------------------------------
%
% Visualize wavepacket in 3 dimensions in momentum (FBR) space 
% using Matlab's builtin isosurface plots
%
%--------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function show_3d_fbr ( obj, state, step )
global expect hamilt info space

if isa(state,'traj')
    prt.error ('Code missing for surface plots of trajectory data in 3D')
end

% Deletes all visible graphics objects from previous time steps
if step>1
    cla 
end

% Loop over wavefunctions/densities
for m=1:hamilt.coupling.n_eqs
    
    % Get densities from trajectories or wavefunctions
    if expect.pop.cha{m}(step)>expect.min_pop
        
        % Get momentum density
        fbr = dvr2fbr(space.dof{1}, state.dvr{m});
        fbr = dvr2fbr(space.dof{2}, fbr);
        fbr = dvr2fbr(space.dof{3}, fbr);
        rho = abs ( fbr ) .^2;
        
        
        % Surface plot of densities
        p = patch(isosurface ( space.fbr{1},space.fbr{2},space.fbr{3}, ...
            rho, obj.scale_fbr/10));
        set(p, ...
            'FaceColor', obj.colors(m,:), ...
            'EdgeColor', 'none')
    end
    
    if m==1
        hold on
    end
    
end
hold off

% Specify view point in terms of azimuth and elevation
view (obj.srf_view(1),obj.srf_view(2))

% Lighting of surfaces
camlight;  
if strcmpi(info.system,'Matlab')
    lighting PHONG;
end
if obj.srf_look(2)
    if strcmpi(info.system,'Matlab')
        lightangle(obj.srf_light(1),obj.srf_light(2));
    elseif strcmpi(info.system,'Octave')
        camlight(obj.srf_light(1),obj.srf_light(2));
end

% Axes, labels, etc
if ~obj.range
    axis ([ ...
        space.dof{1}.fbr_min space.dof{1}.fbr_max ...
        space.dof{2}.fbr_min space.dof{2}.fbr_max ...
        space.dof{3}.fbr_min space.dof{3}.fbr_max ]);
else
    axis ([ ...
        obj.x_min obj.x_max ...
        obj.y_min obj.y_max ...
        obj.z_min obj.z_max ]);
end

set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)

title  ( {info.header1;info.header2} )
xlabel ( [ 'P_{', space.dof{1}.label, '}' ] )
ylabel ( [ 'P_{', space.dof{2}.label, '}' ] )
zlabel ( [ 'P_{', space.dof{3}.label, '}' ] )

end
