%--------------------------------------------------------------------------
% Surface plot of phase space (quasi-)densities and total energy surfaces
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function show_1d_wig ( obj, state, step )
global expect hamilt info space

% Loop over wavefunctions/densities
for m=1:hamilt.coupling.n_eqs
    
    % If population non-negligible
    if expect.pop.cha{m}(step)>expect.min_pop
        
        % Get phase space densities from trajectories or wavefunctions
        if isa(state,'traj')
            x  = space.dof{1}.x_grid;
            y  = space.dof{1}.p_grid;
            n  = space.dof{1}.n_pts;
            dx = space.dof{1}.x_dlt;
            dy = space.dof{1}.p_dlt;
            if strcmpi(info.system,'Matlab') % histcounts2 not available in Octave
                [rho,~] = histcounts2 ( state.pos{1}(state.cha==m), state.mom{1}(state.cha==m), x, y );
            elseif strcmpi(info.system,'Octave') % use of hist3 is discouraged in Matlab
                [rho,~] = hist3 ( [state.pos{1}(state.cha==m), state.mom{1}(state.cha==m)], 'Ctrs',{x y});
            end
            rho = rho / obj.scale_wig;
        elseif isa(state,'wave')
            rho = state.wig{m} / state.wig_max;
        end
    end

    
    % Surface plot of energy functions with|out color coded phase space densities
    if obj.energy
        if expect.pop.cha{m}(step)>expect.min_pop
            
            if isa(state,'traj')
                if strcmpi(info.system,'Matlab')
                    surf ( x(1:n-1)+dx/2, ...
                        y(1:n-1)+dy/2, ...
                        hamilt.tef.wig{m}(1:end-1,1:end-1), ... % not really exact
                        rho' );
                elseif strcmpi(info.system,'Octave')
                    surf ( x(1:n), ...
                        y(1:n), ...
                        hamilt.tef.wig{m}, ... % not really exact
                        rho' );
                end
            elseif isa(state,'wave')
                surf ( ...
                    +space.dvr{1}, ...
                    -space.fbr{1}/2, ...
                    hamilt.tef.wig{m}, ...
                    real(rho) );
            end
        else
            x  = space.dof{1}.x_grid;
            y  = space.dof{1}.p_grid;
            n  = space.dof{1}.n_pts;
            dx = space.dof{1}.x_dlt;
            dy = space.dof{1}.p_dlt;
            surf ( x(1:n-1)+dx/2, ...
                y(1:n-1)+dy/2, ...
                hamilt.tef.wig{m}(1:end-1,1:end-1), ...
                zeros(size(hamilt.tef.wig{m}(1:end-1,1:end-1))) );
        end
        
        % Surface plots of phase space densities
    else
        
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Normalize Wigner distribution: Set maximum to +/- 1
            % Shift individual quasi-densities vertically
            rho = rho + 2*m - 2;
            
            if obj.marginals
                
                if isa(state,'traj') % this includes sub-classes for SHT
                    
                    % Position densities in first row as marginal distribution
                    if strcmpi(info.system,'Matlab') % histcounts not available in Octave
                        [rho_pos,~] = histcounts ( state.pos{1}, x );
                    elseif strcmpi(info.system,'Octave') % use of histc is discouraged in Matlab
                        rho_pos = histc ( state.pos{1}, x );
                    end
                    rho(:,3*n/4) = rho(:,3*n/4) + rho_pos'/obj.scale_dvr;
                    
                    % Momentum densities in first column as marginal distribution
                    if strcmpi(info.system,'Matlab') % histcounts not available in Octave
                        [rho_mom,~] = histcounts ( state.mom{1}, y);
                    elseif strcmpi(info.system,'Octave') % use of histc is discouraged in Matlab
                        rho_mom = histc ( state.mom{1}, y);
                    end
                    rho(1,:) = rho(1,:) + rho_mom/obj.scale_fbr;
                    
                elseif isa(state,'wave')
                    
                    % Position densities in first row as marginal distribution
                    rho_pos = abs ( state.dvr{m} ).^2;
                    rho_pos = min( rho_pos/obj.scale_dvr, 1 );
                    rho(1,:) = rho(1,:) + rho_pos';
                    
                    % Take only inner part of momentum function and double the size
                    state.mom = dvr2fbr(space.dof{1}, state.dvr{m});
                    state.mom = state.mom ( space.dof{1}.n_pts/4+1:3*space.dof{1}.n_pts/4 );
                    state.mom = [state.mom state.mom]';
                    state.mom = reshape ( state.mom, space.dof{1}.n_pts, 1 )/2;
                    
                    % Momentum densities in first column as marginal distribution
                    rho_mom = abs ( state.mom ).^2;
                    rho_mom = min( rho_mom/obj.scale_fbr*4, 1 );
                    rho(:,1) = rho(:,1) + rho_mom;
                    
                end
                
            end
            
            if isa(state,'traj') % this includes sub-classes for SHT
                
                % Dirty trick: Extend z-range from to -1 to 2n-1
                rho(1,1) = 2*hamilt.coupling.n_eqs-1;
                
                % Surface plot of phase space (Wigner) distribution
                surf (  x(1:n-1)+dx/2, y(1:n-1)+dy/2, rho' );
                
            elseif isa(state,'wave')
                
                % Dirty trick: Extend z-range from to -1 to 2n-1
                rho(1,1) = 2*hamilt.coupling.n_eqs-1;
                rho(2,2) = -1;
                
                % Surface plot of phase space (Wigner) distribution
                surf ( space.dvr{1}, -space.fbr{1}/2, rho );
                
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

% Manual versus automatic axis ranges
if ~obj.range
    xyrange = [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
        space.dof{1}.fbr_min/2 space.dof{1}.fbr_max/2 ];
    if obj.energy
        zrange = [ obj.tef_min-obj.tef_delta/10 obj.tef_max+obj.tef_delta/10 ];
    else
        zrange = [ -1 -1+2*hamilt.coupling.n_eqs ];
    end
else
    xyrange = [ obj.x_min obj.x_max obj.y_min obj.y_max ];
    zrange =  [ obj.z_min obj.z_max ];
end
axis ( cat(2, xyrange, zrange) );

% Manual versus automatic color axis
if obj.srf_color
    caxis( [obj.col_min obj.col_max] )
end

% Axis labels
set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)

title  ( {info.header1;info.header2} )
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
ylabel ( [ 'P_{', space.dof{1}.label, '}' ] )

if obj.energy
    if hamilt.coupling.n_eqs==1
        zlabel ( 'E(R,P)' )
    else
        if strcmpi ( hamilt.coupling.represent,'adi')
            zlabel ( 'E_{adi}(R,P)' )
        elseif strcmpi ( hamilt.coupling.represent,'dia')
            zlabel ( 'E_{dia}(R,P)' )
        end
    end
else
    if hamilt.coupling.n_eqs==1
        zlabel ( '\rho(R,P)' )
    else
        if strcmpi ( hamilt.coupling.represent,'adi')
            zlabel ( '\rho_{adi}(R,P)' )
        elseif strcmpi ( hamilt.coupling.represent,'dia')
            zlabel ( '\rho_{dia}(R,P)' )
        end
    end
end

end



