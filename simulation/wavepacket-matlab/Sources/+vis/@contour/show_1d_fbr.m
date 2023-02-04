%--------------------------------------------------------------------------
% Plot momentum density and kinetic energy curve (from wavefunctions)
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function show_1d_fbr ( obj, state, step )
global expect info hamilt space

%% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs
    
    % If population not too small
    if expect.pop.cha{m}(step) > expect.min_pop
        
        % Horizontal offset (= kinetic + potential energy)
        if obj.energy
            offset = expect.pot.cha{m}(step) + expect.kin.cha{m}(step);
            scaling = obj.kin_delta;
            y_min = obj.kin_min-obj.kin_delta/10;
            y_max = obj.kin_max+obj.kin_delta/10;
            
        else
            offset = 0;
            scaling = 1;
            y_min = -0.1;
            y_max = +1.1;
        end
        
        
        if isa(state,'traj') % this includes sub-classes for SHT
            
            y  = space.dof{1}.p_grid;
            n  = space.dof{1}.n_pts;
            dy = space.dof{1}.p_dlt;
            
            if strcmpi(info.system,'Matlab') % histcounts not available in Octave
                [rho,~] = histcounts (state.mom{1}(state.cha==m), y);
                h = plot ( rho*scaling/obj.scale_fbr+offset, y(1:n-1)+dy/2 );
            elseif strcmpi(info.system,'Octave') % use of histc is discouraged in Matlab
                rho = histc (state.mom{1}(state.cha==m), y);
                h = plot ( rho*scaling/obj.scale_fbr+offset, y(1:n) );
            end
            
            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     obj.colors(m,:), ...
                'LineWidth', obj.l_thick)
            
        elseif isa(state,'wave')
            
            state.fbr{m} = dvr2fbr(space.dof{1}, state.dvr{m});
            
            % Get density from wavefunction
            rho = abs ( state.fbr{m}(space.dof{1}.n_pts/4+1 : 3*space.dof{1}.n_pts/4) ) .^2;
            
            % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
            phi = angle ( state.fbr{m}(space.dof{1}.n_pts/4+1 : 3*space.dof{1}.n_pts/4) ) / (2*(pi+0.001)) + 1/2;
            
            if strcmpi(info.system,'Matlab')
                vis.styles.plot_color ( space.fbr{1}(space.dof{1}.n_pts/4+1 : 3*space.dof{1}.n_pts/4), ...
                    rho*scaling/obj.scale_fbr, ...
                    phi, ...
                    obj.colors(m,:), ...
                    obj.l_thick, ...
                    offset, ...
                    1 )
            elseif strcmpi(info.system,'Octave')
                h = plot ( rho*scaling/obj.scale_fbr + offset, space.fbr{1}(space.dof{1}.n_pts/4+1 : 3*space.dof{1}.n_pts/4) );
                set(h, ...
                    'LineStyle', obj.patterns{1}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thick )
            end
            
        end
        
        % Plot vertical base line
        if obj.energy
            h = line ( [1 1]*offset, [space.dof{1}.fbr_min space.dof{1}.fbr_max] );
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

    h = plot ( space.dof{1}.kin, plusminus * space.fbr{1} );
    set(h, ...
        'LineStyle', obj.patterns{1}, ...
        'Color',     'black', ...
        'LineWidth', obj.l_thin)
end

%% Axes, labels, etc
if ~obj.range
    axis ( [ y_min y_max space.dof{1}.fbr_min/2 space.dof{1}.fbr_max/2 ] )
else
    axis ( [ y_min y_max obj.y_min obj.y_max ] )
end

set(gca, ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', ...
    'LineWidth',     obj.l_thick, ...
    'FontName',      obj.f_name, ...
    'FontSize',      obj.f_large, ...
    'FontWeight',    obj.f_heavy)

if obj.energy
    xlabel ( 'T(P)' )
else
    
    if hamilt.coupling.n_eqs==1
        xlabel ( '\rho(P)' )
    else
        if strcmpi ( hamilt.coupling.represent,'adi' )
            xlabel ( '\rho_{adi}(P)' )
        elseif strcmpi ( hamilt.coupling.represent,'dia' )
            xlabel ( '\rho_{dia}(P)' )
        end
    end
    
end
ylabel ( [ 'P_{', space.dof{1}.label, '}' ] )

end

