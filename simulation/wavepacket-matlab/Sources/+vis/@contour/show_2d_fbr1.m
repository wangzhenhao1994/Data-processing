%------------------------------------------------------------
% Densities projected on first coordinate (vertical plot)
% Momentum space
%------------------------------------------------------------
function show_2d_fbr1 ( obj, state, step )      
global expect hamilt info space
persistent rho_max

% Preallocation
rho = cell (hamilt.coupling.n_eqs,1);

% Get reduced densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step)>expect.min_pop
        switch(lower(class(state)))
            case ('traj')
                if strcmpi(info.system,'Matlab') % histcounts not available in Octave
                    [c,~] = histcounts (state.mom{1}(state.cha==m), space.dof{1}.p_grid);
                elseif strcmpi(info.system,'Octave') % use of histc is discouraged in Matlab
                    c = histc (state.mom{1}(state.cha==m), space.dof{1}.p_grid);
                end
                rho{m} = c';
            case ('wave')
                fbr = dvr2fbr(space.dof{1}, state.dvr{m});
                fbr = dvr2fbr(space.dof{2}, fbr);
                rho{m} = sum( abs(fbr).^2, 2 );
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
            y  = space.dof{1}.p_grid;
            n  = space.dof{1}.n_pts;
            dy = space.dof{1}.p_dlt;
            if strcmpi(info.system,'Matlab')
                h = plot ( y(1:n-1)+dy/2, rho{m}/rho_max );
            elseif strcmpi(info.system,'Octave')
                h = plot ( y(1:n), rho{m}/rho_max );
            end
        elseif isa(state,'wave')
            h = plot ( space.dof{1}.p_grid,rho{m}/rho_max );
        end
        set(h, ...
            'LineStyle', obj.patterns{1}, ...
            'LineWidth', obj.l_thick, ...
            'Color',     obj.colors(m,:))
    end
end

if ~obj.range
    axis ( [ space.dof{1}.fbr_min space.dof{1}.fbr_max 0 1 ] )
else
    axis ( [ obj.x_min obj.x_max 0 1 ] )
end

% Axes and labels
set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)


xlabel ( [ 'P_{', space.dof{1}.label, '}' ] )
ylabel ( [ '\rho (P_{', space.dof{1}.label, '})' ] )


