%--------------------------------------------------------------------------
%
% Calculate and visualize flux density (1D)
%
% J = -i*hbar/2 * [       psi^* * grad (psi) - psi * grad (psi^*) ]
% 
%   =    hbar * Im [      psi^* * grad (psi) ]
%
%   =    hbar * Re [ -i * psi^* * grad (psi) ]
%
%   =           Re [      psi^* *   P  (psi) ]
%
% see Cohen-Tannoudji, Chapter III, Complement D, Eq. D-17 (p.239)
%
% normalization with 1/mass ?!?
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function show_2d ( obj, state, step )

global expect hamilt info space

% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs
    
    % Plot dotted contours of total energy functions
    if obj.energy
        if ~isempty(hamilt.pot{m,m}.dvr)
            [~,h] = contour ( space.dvr{1} , ...
                space.dvr{2}, ...
                hamilt.pot{m,m}.dvr);
        else
            [~,h] = contour ( space.dvr{1} , ...
                space.dvr{2}, ...
                zeros(size(space.dvr{1})));
        end
        set(h, ...
            'LevelList', linspace( obj.pot_min, obj.pot_max, obj.cnt_nlev(2) ), ...
            'LineStyle', obj.patterns{3}, ...
            'LineColor', obj.colors(m,:), ...
            'LineWidth', obj.l_thin)
    end
    
    % If populations not too small, plot solid contours of densities
    if expect.pop.cha{m}(step)>expect.min_pop
        
        % Contour levels
        if obj.cnt_levels
            cmin = obj.cnt_min;
            cmax = obj.cnt_max;
        else
            cmin = 1 / obj.cnt_nlev(1);
            cmax = 1;
        end
        
        if isa (state,'traj')
            
            % Contour plot of position densities
            x  = space.dof{1}.x_grid;
            y  = space.dof{2}.x_grid;
            nx = space.dof{1}.n_pts;
            ny = space.dof{2}.n_pts;
            dx = space.dof{1}.x_dlt;
            dy = space.dof{2}.x_dlt;
            if strcmpi(info.system,'Matlab') % histcounts2 not available in Octave
                [rho,~] = histcounts2 ( state.pos{1}(state.cha==m), state.pos{2}(state.cha==m), x, y );
                [~,h] = contour ( x(1:nx-1)+dx/2, y(1:ny-1)+dy/2, rho'/obj.scale_dvr );
            elseif strcmpi(info.system,'Octave') % use of hist3 is discouraged in Matlab
                [rho,~] = hist3 ( [state.pos{1}(state.cha==m), state.pos{2}(state.cha==m)], 'Ctrs',{x y} );
                [~,h] = contour ( x(1:nx), y(1:ny), rho'/obj.scale_dvr );
            end
            set(h, ...
                'LevelList', linspace(cmin, cmax, 5), ...
                'LineStyle', obj.patterns{2}, ...
                'LineColor', obj.colors(m,:), ...
                'LineWidth', obj.l_thin)
            
            % Quiver plot of momentum vectors
            hold on
            quiver ( ...
                state.pos{1}(state.cha==m), ...
                state.pos{2}(state.cha==m), ...
                state.mom{1}(state.cha==m), ...
                state.mom{2}(state.cha==m) );
            hold off
            
        elseif isa (state,'wave')
            
            % Contour plot of position densities
            [~,h] = contour ( ...
                space.dvr{1},    ...
                space.dvr{2}, ...
                abs(state.dvr{m}).^2 / obj.scale_dvr);
            set(h, ...
                'LevelList', linspace(cmin, cmax, 5), ...
                'LineStyle', obj.patterns{2}, ...
                'LineColor', obj.colors(m,:), ...
                'LineWidth', obj.l_thin)
            
            % Quiver plot of flux density
            hold on
            psi_mom = momentum(space.dof{1}, state.dvr{1});
            j1 = real ( conj(state.dvr{m}) .* psi_mom );
            psi_mom = momentum(space.dof{2}, state.dvr{1});
            j2 = real ( conj(state.dvr{m}) .* psi_mom );
            quiver ( ...
                space.dvr{1}, ...
                space.dvr{2}, ...
                j1, ...
                j2 );
            hold off
            
        end
    end
    
end

if ~obj.range
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
            space.dof{2}.dvr_min space.dof{2}.dvr_max ] )
else
    axis( [ obj.x_min obj.x_max obj.y_min obj.y_max ] )
end

set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)

title ( {info.header1;info.header2} )
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )

