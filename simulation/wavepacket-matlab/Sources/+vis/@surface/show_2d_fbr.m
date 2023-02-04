%--------------------------------------------------------------
% Surface plot of momentum (FBR) densities in two dimensions
%--------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function show_2d_fbr ( obj, state, step )
global expect hamilt info space

% Prepare histograms for trajectory data
if isa(state,'traj')
    x  = space.dof{1}.p_grid;
    y  = space.dof{2}.p_grid;
    nx = space.dof{1}.n_pts;
    ny = space.dof{2}.n_pts;
    dx = space.dof{1}.p_dlt;
    dy = space.dof{2}.p_dlt;
end

% Loop over wavefunctions/densities
for m=1:hamilt.coupling.n_eqs
    
    % If population non-negligible
    if expect.pop.cha{m}(step)>expect.min_pop
        
        % Get position densities from trajectories or wavefunctions
        if isa(state,'traj')
            
            if strcmpi(info.system,'Matlab') % histcounts2 not available in Octave
                [rho,~] = histcounts2 ( state.mom{1}(state.cha==m), state.mom{2}(state.cha==m), x, y );
            elseif strcmpi(info.system,'Octave') % use of hist3 is discouraged in Matlab
                [rho,~] = hist3 ( [state.mom{1}(state.cha==m), state.mom{2}(state.cha==m)], 'Ctrs',{x y} );
            end
            
        elseif isa(state,'wave')
            
            fbr = dvr2fbr(space.dof{1}, state.dvr{m});
            fbr = dvr2fbr(space.dof{2}, fbr);
            rho = abs(fbr).^2;
        end
        rho = rho / obj.scale_fbr;
    end
    
    % Surface plots of potentials with|out color coded densities
    if obj.energy && isa(space.dof{1}, 'dof.fft') ...
            && isa(space.dof{2}, 'dof.fft')
        if expect.pop.cha{m}(step)>expect.min_pop
            
            if isa(state,'traj')
                kinetic = space.dof{1}.kin + space.dof{2}.kin;
                surf ( ...
                    x(1:nx-1)+dx/2, ...
                    y(1:ny-1)+dy/2, ...
                    kinetic(1:end-1,1:end-1)', ... % not really exact
                    rho' );
                
            elseif isa(state,'wave')
                
                surf ( ...
                    space.fbr{1}, ...
                    space.fbr{2}, ...
                    space.dof{1}.kin + space.dof{2}.kin, ...
                    rho );
            end
        else
            
            if m==1
                surf ( space.fbr{1}, ...
                       space.fbr{2}, ...
                       space.dof{1}.kin + space.dof{2}.kin );
            end
            
        end
        
    % Surface plots of densities only (and color-coded phases)
    else
        
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Normalized densities: Set maximum to + 1
            % Shift individual densities vertically
            rho = rho + 2*m - 2;
            
            % Surface plot of momentum densities
            if isa(state,'traj')
                
                if strcmpi(info.system,'Matlab')
                    surf ( x(1:nx-1)+dx/2, y(1:ny-1)+dy/2, rho' );
                elseif strcmpi(info.system,'Octave')
                    surf ( x(1:nx), y(1:ny), rho' );
                end
                
            elseif isa(state,'wave')
                
                % Dirty trick: Extend z-range from to -1 to 2*n-1
                rho( 1,  1 ) = 2*hamilt.coupling.n_eqs-1;
                rho(end,end) = -1;
                
                % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
                phi = angle ( fbr ) / (2*(pi+0.001)) + 1/2;
                phi(rho<0.01) = 0;
                
                % Surface plot of densities, with phases as colors
                surf ( space.fbr{1}, space.fbr{2}, rho, phi );
                
            end
            
        end
        
    if m==1
        hold on
    end

end
hold off

% Specify view point in terms of azimuth and elevation
view (obj.srf_view(1),obj.srf_view(2))

% Select color map
colormap(obj.col_map);

% Shading of surface(s)
if obj.srf_look(1)
    shading interp;
end

% Lighting of surfaces
if obj.srf_look(2)
    if strcmpi(info.system,'Matlab')
        lightangle(obj.srf_light(1),obj.srf_light(2));
    elseif strcmpi(info.system,'Octave')
        camlight(obj.srf_light(1),obj.srf_light(2));
end

% Axes, labels, etc
set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)

if ~obj.range
    xyrange = [ space.dof{1}.fbr_min space.dof{1}.fbr_max ...
        space.dof{2}.fbr_min space.dof{2}.fbr_max ];
    if obj.energy
        zrange = [ obj.kin_min-obj.kin_delta/10 obj.kin_max+obj.kin_delta/10 ];
    else
        zrange = [ -1 -1+2*hamilt.coupling.n_eqs ];
    end
else
    xyrange = [ obj.x_min obj.x_max obj.y_min obj.y_max ];
    zrange  = [ obj.z_min obj.z_max ];
end
axis ( cat(2, xyrange, zrange) );

% Manual versus automatic color axis
if obj.srf_color
    caxis( [obj.col_min obj.col_max] )
end

title  ( {info.header1;info.header2} )
xlabel ( [ 'P_{', space.dof{1}.label, '}' ] )
ylabel ( [ 'P_{', space.dof{2}.label, '}' ] )

if obj.energy && isa(space.dof{1}, 'dof.fft') ...
        && isa(space.dof{2}, 'dof.fft')
    zlabel ( 'T(P)' )
else
    if hamilt.coupling.n_eqs==1
        zlabel ( '\rho(P) ' )
    else
        if strcmpi ( hamilt.coupling.represent,'adi' )
            zlabel ( '\rho_{adi}(P)' )
        elseif strcmpi ( hamilt.coupling.represent,'dia' )
            zlabel ( '\rho_{dia}(P)' )
        end
    end
end

end


