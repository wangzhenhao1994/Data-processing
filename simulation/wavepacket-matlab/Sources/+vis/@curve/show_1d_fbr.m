%----------------------------------------------------------
% Plot momentum density and kinetic energy curve
%----------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2009 Ulf Lorenz
%
% see the README file for license details.

function show_1d_fbr ( obj, state, step )
global expect hamilt info space

%% Loop over (coupled) channels
for  m=1:hamilt.coupling.n_eqs
    
    % If population not too small
    if ( expect.pop.cha{m}(step) > expect.min_pop )
        
        % Horizontal offset (= kinetic + potential energy)
        if obj.energy
            offset = expect.kin.cha{m}(step) + expect.pot.cha{m}(step);
            scaling = obj.kin_delta;
            y_min = obj.kin_min-obj.kin_delta/10;
            y_max = obj.kin_max+obj.kin_delta/10;
            
        else
            offset = 0;
            scaling = 1;
            y_min = -0.1;
            y_max = +1.1;
            
        end
        
        
        if isa(state,'traj')  % this includes sub-classes for SHT
            
            y  = space.dof{1}.p_grid;
            n  = space.dof{1}.n_pts;
            dy = space.dof{1}.p_dlt;
            
            if strcmpi(info.system,'Matlab') % histcounts not available in Octave
                [rho,~] = histcounts (state.mom{1}(state.cha==m), y);
                h = plot ( y(1:n-1)+dy/2, rho*scaling/obj.scale_fbr+offset);
            elseif strcmpi(info.system,'Octave') % use of histc is discouraged in Matlab
                rho = histc (state.mom{1}(state.cha==m), y);
                h = plot ( y(1:n), rho*scaling/obj.scale_fbr+offset);
            end
            
            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     obj.colors(m,:), ...
                'LineWidth', obj.l_thick)
            
        elseif isa (state,'wave')
            
            state.fbr{m} = dvr2fbr(space.dof{1}, state.dvr{m});
            
            % Get density from wavefunction
            rho = abs ( state.fbr{m} ) .^2;
            
            % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
            phi = angle ( state.fbr{m} ) / (2*(pi+0.001)) + 1/2;
            
            if strcmpi(info.system,'Matlab')
                vis.styles.plot_color ( space.fbr{1}, ...
                    rho*scaling/obj.scale_fbr, ...
                    phi, ...
                    obj.colors(m,:), ...
                    obj.l_extra, ...
                    offset, ...
                    0 )
            elseif  strcmpi(info.system,'Octave')
                h = plot ( space.fbr{1}, rho*scaling/obj.scale_fbr + offset);
                set(h, ...
                    'LineStyle', obj.patterns{1}, ...
                    'LineWidth', obj.l_thick, ...
                    'Color',     obj.colors(m,:) )
            end
            
        end
        
        % Plot horizontal base line
        if obj.energy
            h = line ( [space.dof{1}.fbr_min space.dof{1}.fbr_max], [1 1]*offset );
            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     obj.colors(m,:), ...
                'LineWidth', obj.l_thin)
        end
        
    end
    
end

%% Plot kinetic energy curve
if obj.energy
    if isa(space.dof{1}, 'dof.fft')
        plusminus = -1;
    else
        plusminus = +1;
    end
    h = plot ( plusminus * space.fbr{1}, space.dof{1}.kin);
    set(h, ...
        'LineStyle', obj.patterns{1}, ...
        'Color',     'black', ...
        'LineWidth', obj.l_thin)
end

%% Axes, labels, etc
if ~obj.range
    axis ( [ space.dof{1}.fbr_min space.dof{1}.fbr_max y_min y_max ] )
else
    if ~isempty(obj.y_min) && ~isempty(obj.y_max) 
        axis ( [ obj.x_min obj.x_max obj.y_min obj.y_max ] )
    else
        axis ( [ obj.x_min obj.x_max y_min y_max ] )
    end
end

set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)

xlabel ( [ 'P_{', space.dof{1}.label, '}' ] )
title ( {info.header1;info.header2} )

if obj.energy
    ylabel ( 'T(P)' )
else
    if hamilt.coupling.n_eqs==1
        ylabel ( '\rho(P)' )
    else
        if strcmpi ( hamilt.coupling.represent,'adi' )
            ylabel ( '\rho_{adi}(P)' )
        elseif strcmpi ( hamilt.coupling.represent,'dia' )
            ylabel ( '\rho_{dia}(P)' )
        end
    end
    
end

end
 
