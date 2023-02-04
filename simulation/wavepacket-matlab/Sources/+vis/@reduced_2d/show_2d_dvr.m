%-------------------------------------------------------------
% Contour plots of reduced density matrices for (2*k-1)-th and 
% (2*k)-th degree of freedom in position (R_k) representation
%-------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20xy Burkhard Schmidt's group
%
% see the README file for license details.

function show_2d_dvr ( obj, state, step, k )       
global info space expect hamilt

% Initialize calculations of densities and their maxima
rho = cell(hamilt.coupling.n_eqs,1);
if step==1
    obj.rho_max(k) = 0;
end

% Main switch: Q/M versus Q/C
if isa(state,'traj') % this includes sub-classes for SHT
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Get 2-dim histograms from trajectory bundles
            x  = space.dof{2*k-1}.x_grid;
            y  = space.dof{2*k  }.x_grid;
            nx = space.dof{2*k-1}.n_pts;
            ny = space.dof{2*k  }.n_pts;
            dx = space.dof{2*k-1}.x_dlt;
            dy = space.dof{2*k  }.x_dlt;
            
            if strcmpi(info.system,'Matlab') % histcounts2 not available in Octave
                [rho{m},~] = histcounts2 ( state.pos{2*k-1}(state.cha==m), state.pos{2*k}(state.cha==m), x, y );
            elseif strcmpi(info.system,'Octave') % use of hist3 is discouraged in Matlab
                [rho{m},~] = hist3 ( [state.pos{2*k-1}(state.cha==m), state.pos{2*k}(state.cha==m)], 'Ctrs',{x y} );
            end
            
            % Find maximal density at first step
            if step==1
                obj.rho_max(k) = max ( obj.rho_max(k), max(max(rho{m})) );
            end
            
        end
    end
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Create curve plots
            if strcmpi(info.system,'Matlab')
                [~,h] = contour ( x(1:nx-1)+dx/2, y(1:ny-1)+dy/2, rho{m}'/obj.rho_max(k) );
            elseif strcmpi(info.system,'Octave')
                [~,h] = contour ( x(1:nx), y(1:ny), rho{m}'/obj.rho_max(k) );
            end
            
            set(h, ...
                'LevelList', linspace(0, 1, obj.cnt_nlev(1)), ...   % use even number of contours to avoid zero!
                'LineStyle', obj.patterns{1}, ...
                'LineColor', obj.colors(m,:), ...
                'LineWidth', obj.l_thin)
        end
    end
    
elseif isa(state,'wave')
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Find maximal density at first step
            if step==1
                obj.rho_max(k) = max ( obj.rho_max(k), max(max(state.redu{m,k})) );
            end
            
        end
    end
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Create contour plots
            [~,h] = contour ( ...
                space.dof{2*k-1}.x_grid, ...
                space.dof{2*k  }.x_grid, ...
                state.redu{m,k});
            set(h, ...
                'LevelList', linspace(0, obj.rho_max(k), obj.cnt_nlev), ...   % use even number of contours to avoid zero!
                'LineStyle', obj.patterns{1}, ...
                'LineColor', obj.colors(m,:), ...
                'LineWidth', obj.l_thick)
        end
    end
    
end

% Optionally setting plot ranges "manually"
if ~obj.range
    axis ( [ space.dof{2*k-1}.dvr_min space.dof{2*k-1}.dvr_max ...
        space.dof{2*k}.dvr_min space.dof{2*k}.dvr_max ] )
else
    axis ( [ obj.x_min obj.x_max obj.y_min obj.y_max ] )
end

% Thick lines and heavy fonts
set(gca, ...
	'LineWidth',  obj.l_thick, ...
	'FontName',   obj.f_name, ...
	'FontSize',   obj.f_large, ...
	'FontWeight', obj.f_heavy)

% Place the header approximately in the middle
if k==2
    title ( {info.header1;info.header2} )
end

% Labels of x,y axes
xlabel ( ['R_{', space.dof{2*k-1}.label, '}'] )
ylabel ( ['R_{', space.dof{2*k  }.label, '}'] )


% Negative imaginary potential (as absorbing boundary conditions)
if isa(state,'wave') && isfield (hamilt,'nip')
    for  m=1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.nip{m}.dvr)
            
            % Left border
            if hamilt.nip{m}.min(2*k-1) > space.dof{2*k-1}.dvr_min 
                h = line ( [ hamilt.nip{m}.min(2*k-1) hamilt.nip{m}.min(2*k-1) ], ...
                    [ space.dof{2*k}.dvr_min space.dof{2*k}.dvr_max ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
            
            % Lower border
            if hamilt.nip{m}.min(2*k) > space.dof{2*k}.dvr_min 
                h = line ( [ space.dof{2*k-1}.dvr_min space.dof{2*k-1}.dvr_max ], ...
                    [ hamilt.nip{m}.min(2*k) hamilt.nip{m}.min(2*k) ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
            
            % Right border
            if hamilt.nip{m}.max(2*k-1) < space.dof{2*k-1}.dvr_max
                h = line ( [ hamilt.nip{m}.max(2*k-1) hamilt.nip{m}.max(2*k-1) ], ...
                    [ space.dof{2*k}.dvr_min space.dof{2*k}.dvr_max ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
            
            % Upper border
            if hamilt.nip{m}.max(2*k) < space.dof{2*k}.dvr_max
                h = line ( [ space.dof{2*k-1}.dvr_min space.dof{2*k-1}.dvr_max ], ...
                    [ hamilt.nip{m}.max(2*k) hamilt.nip{m}.max(2*k) ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
        end
    end
end
