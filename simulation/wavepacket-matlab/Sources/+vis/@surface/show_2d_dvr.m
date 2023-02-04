%--------------------------------------------------------------------------
% Surface plot of position (DVR) densities in two dimensions
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2019 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function show_2d_dvr ( obj, state, step )
global expect hamilt info space

% Prepare histograms for trajectory data
if isa(state,'traj')
    x  = space.dof{1}.x_grid;
    y  = space.dof{2}.x_grid;
    nx = space.dof{1}.n_pts;
    ny = space.dof{2}.n_pts;
    dx = space.dof{1}.x_dlt;
    dy = space.dof{2}.x_dlt;
end

% Loop over wavefunctions/densities
for m=1:hamilt.coupling.n_eqs
    
    % Get densities from trajectories or wavefunctions
    if expect.pop.cha{m}(step)>expect.min_pop
        
        if isa(state,'traj')
            
            if strcmpi(info.system,'Matlab') % histcounts2 not available in Octave
                [rho,~] = histcounts2 ( state.pos{1}(state.cha==m), state.pos{2}(state.cha==m), x, y );
            elseif strcmpi(info.system,'Octave') % use of hist3 is discouraged in Matlab
                [rho,~] = hist3 ( [state.pos{1}(state.cha==m), state.pos{2}(state.cha==m)], 'Ctrs',{x y} );
            end
         
        elseif isa(state,'wave')
            
            rho = abs(state.dvr{m}).^2;
            
        end
        rho = rho / obj.scale_dvr;
    end
    
    % Plots of potential energy surfaces
    if obj.energy
        
        % Surfaces with densities color coded
        if expect.pop.cha{m}(step)>expect.min_pop
            
            if isa(state,'traj')
                if strcmpi(info.system,'Matlab')
                    surf (  ...
                        x(1:nx-1)+dx/2, ...
                        y(1:ny-1)+dy/2, ...
                        hamilt.pot{m,m}.dvr(1:end-1,1:end-1)',... % not really exact
                        rho' );
                elseif strcmpi(info.system,'Octave')
                    surf (  ...
                        x(1:nx), ...
                        y(1:ny), ...
                        hamilt.pot{m,m}.dvr(1:end,1:end)',... % not really exact
                        rho' );
                end
            elseif isa(state,'wave')
                if ~isempty(hamilt.pot{m,m}.dvr)
                    surf ( ...
                        space.dvr{1}, ...
                        space.dvr{2}, ...
                        hamilt.pot{m,m}.dvr, ...
                        rho );
                else
                    surf ( ...
                        space.dvr{1}, ...
                        space.dvr{2}, ...
                        zeros(size(space.dvr{1})), ...
                        rho );
                end
            end
            
            % Surfaces only  if (practically) no density
        else
            if isa(state,'traj')
                surf (  ...
                    x(1:nx-1)+dx/2, ...
                    y(1:ny-1)+dy/2, ...
                    hamilt.pot{m,m}.dvr(1:end-1,1:end-1)',...
                    zeros (nx-1,ny-1)') % not really exact
                
            elseif isa(state,'wave')
                surf ( ...
                    space.dvr{1}, ...
                    space.dvr{2}, ...
                    hamilt.pot{m,m}.dvr, ...
                    zeros (space.dof{1}.n_pts,space.dof{2}.n_pts));
                
            end
        end
        
        % Surface plots of densities only (and color-coded phases)
    else
        
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Normalized densities: Set maximum to + 1
            % Shift individual densities vertically
            rho = rho + 2*m - 2;
            
            % Dirty trick: Extend z-range from to -1 to 2*n-1
            rho( 1,  1) = 2*hamilt.coupling.n_eqs-1;
            rho(end,end) = -1;
            
            % Surface plot of position densities
            if isa(state,'traj')
                
                if strcmpi(info.system,'Matlab')
                    surf (  x(1:nx-1)+dx/2, y(1:ny-1)+dy/2, rho' );
                elseif strcmpi(info.system,'Octave')
                    surf (  x(1:nx), y(1:ny), rho' );
                end
                
            elseif isa(state,'wave')
                
                
                % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
                phi = angle ( state.dvr{m} ) / (2*(pi+0.001)) + 1/2;
                phi(rho<0.01) = 0;
                
                % Surface plot of densities, with phases as colors
                surf ( space.dvr{1}, space.dvr{2}, rho, phi );
                
            end
            
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
    xyrange = [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
        space.dof{2}.dvr_min space.dof{2}.dvr_max ];
    if obj.energy
        zrange = [ obj.pot_min-obj.pot_delta/10 obj.pot_max+obj.pot_delta/10 ];
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
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )

if obj.energy
    
    zlabel ( 'Potential energy' )
    if hamilt.coupling.n_eqs==1
        zlabel ( 'V(R)' )
    else
        if strcmpi ( hamilt.coupling.represent,'adi' )
            zlabel ( 'V_{adi}(R)' )
        elseif strcmpi ( hamilt.coupling.represent,'dia' )
            zlabel ( 'V_{dia}(R)' )
        end
    end
else
    
    if hamilt.coupling.n_eqs==1
        zlabel ( '\rho(R) ' )
    else
        if strcmpi ( hamilt.coupling.represent,'adi' )
            zlabel ( '\rho_{adi}(R)' )
        elseif strcmpi ( hamilt.coupling.represent,'dia' )
            zlabel ( '\rho_{dia}(R)' )
        end
    end
end

end

