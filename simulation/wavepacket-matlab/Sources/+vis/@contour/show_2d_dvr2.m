%------------------------------------------------------------
% Densities projected on second coordinate (horizontal plot)
%------------------------------------------------------------
function show_2d_dvr2 ( obj, state, step )       
global expect hamilt info space
persistent rho_max

% Preallocate
rho  = cell (hamilt.coupling.n_eqs,1);

% Get reduced densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step)>expect.min_pop
        if isa(state,'traj')
            if strcmpi(info.system,'Matlab') % histcounts not available in Octave
                [c,~] = histcounts (state.pos{2}(state.cha==m), space.dof{2}.x_grid);
            elseif strcmpi(info.system,'Octave') % use of histc is discouraged in Matlab
                c = histc (state.pos{2}(state.cha==m), space.dof{2}.x_grid);
            end
            rho{m} = c';
        elseif isa(state,'wave')
            rho{m} = sum( abs(state.dvr{m}).^2 .* space.weight, 1 );
        end
    end
end
        
% Get maximum density
% if step==1
rho_max = 0;
for  m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step)>expect.min_pop
        rho_max = max ( rho_max, max(rho{m}) );
    end
end
% end

% Plot projected densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step)>expect.min_pop
        if isa(state,'traj')
            x  = space.dof{2}.x_grid;
            n  = space.dof{2}.n_pts;
            dx = space.dof{2}.x_dlt;
            if strcmpi(info.system,'Matlab')
                h = plot ( rho{m}/rho_max, x(1:n-1)+dx/2 );
            elseif strcmpi(info.system,'Octave')
                h = plot ( rho{m}/rho_max, x(1:n) );
            end
        elseif isa(state,'wave')
            h = plot ( rho{m}/rho_max, space.dof{2}.x_grid );
        end
        
        set(h, ...
            'LineStyle', obj.patterns{1}, ...
            'LineWidth', obj.l_thick, ...
            'Color',     obj.colors(m,:))

    end
end

% Axes and labels
set(gca, ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', ...
    'LineWidth',     obj.l_thick, ...
    'FontName',      obj.f_name, ...
    'FontSize',      obj.f_large, ...
    'FontWeight',    obj.f_heavy)

xlabel ( [ '\rho (R_{', space.dof{2}.label, '})' ] )
ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )
if ~obj.range
    axis ( [ 0 1 space.dof{2}.dvr_min space.dof{2}.dvr_max ] )
else
    axis ( [ 0 1 obj.y_min obj.y_max ] )
end

% Absorbing boundary conditions
if isa(state,'wave') && isfield (hamilt,'nip')
    for  m=1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.nip{m}.dvr)
            if hamilt.nip{m}.min(2) > space.dof{2}.dvr_min % Lower
                h = line ( [ 0 1 ], [ hamilt.nip{m}.min(2) hamilt.nip{m}.min(2) ] );

                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'LineWidth', obj.l_thin, ...
                    'Color',     obj.colors(m,:))

            end
            
            if hamilt.nip{m}.max(2) < space.dof{2}.dvr_max % Upper
                h = line ( [ 0 1 ], [ hamilt.nip{m}.max(2) hamilt.nip{m}.max(2) ] );

                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'LineWidth', obj.l_thin, ...
                    'Color',     obj.colors(m,:))

            end
        end
    end
end

