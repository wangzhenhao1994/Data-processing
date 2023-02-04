%------------------------------------------------------------
% Plot 2-dim momentum densities and kinetic energy surfaces
%------------------------------------------------------------
function show_2d_fbr ( obj, state, step )       
global expect hamilt info space

% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs
    
    % Plot "mean trajectory" as solid curve in position space
    % Expectation value of R_1 vs. R_2
    if obj.expect
        plots.mask.ind{m} = find ( expect.pop.cha{m} > expect.min_pop );
        h = plot ( expect.mom{1}.cha{m}(plots.mask.ind{m}), ...
             expect.mom{2}.cha{m}(plots.mask.ind{m}));
        set(h, ...
            'LineStyle', obj.patterns{1}, ...
            'LineWidth', obj.l_thick, ...
            'Color',     obj.colors(m,:))
    end
       
    % Plot dotted contours of kinetic enery functions, for non-FFt grids only
    if obj.energy && isa(space.dof{1}, 'dof.fft') ...
            && isa(space.dof{2}, 'dof.fft')
        [~,h] = contour ( space.fbr{1} , ...
            space.fbr{2}, ...
            space.dof{1}.kin + space.dof{2}.kin);
        set(h, ...
            'LevelList', linspace( obj.kin_min, obj.kin_max, obj.cnt_nlev(2) ), ...
            'LineStyle', obj.patterns{3}, ...
            'LineWidth', obj.l_thin, ...
            'LineColor', obj.colors(m,:))

    end
    
    % If populations not too small, plot solid contours of densities
    if expect.pop.cha{m}(step)>expect.min_pop
        
        if obj.cnt_levels
            cmin = obj.cnt_min;
            cmax = obj.cnt_max;
        else
            cmin = 1 / obj.cnt_nlev(1);
            cmax = 1;
        end
        
        if isa(state,'traj')
            x  = space.dof{1}.p_grid;
            y  = space.dof{2}.p_grid;
            nx = space.dof{1}.n_pts;
            ny = space.dof{2}.n_pts;
            dx = space.dof{1}.p_dlt;
            dy = space.dof{2}.p_dlt;
            
            if strcmpi(info.system,'Matlab') % histcounts2 not available in Octave
                [rho,~] = histcounts2 ( state.mom{1}(state.cha==m), state.mom{2}(state.cha==m), x, y );
                [~,h] = contour (  x(1:nx-1)+dx/2, y(1:ny-1)+dy/2, rho'/obj.scale_dvr );
            elseif strcmpi(info.system,'Octave') % use of hist3 is discouraged in Matlab
                [rho,~] = hist3 ( [state.mom{1}(state.cha==m), state.mom{2}(state.cha==m)], 'Ctrs',{x y} );
                [~,h] = contour (  x(1:nx), y(1:ny), rho'/obj.scale_dvr );
            end
            
        elseif isa(state,'wave')
            fbr = dvr2fbr(space.dof{1}, state.dvr{m});
            fbr = dvr2fbr(space.dof{2}, fbr);
            
            switch lower(obj.complex)
                
                case 'abs2'
                    [~,h] = contour ( ...
                        space.fbr{1}, ...
                        space.fbr{2}, ...
                        abs(fbr).^2 / obj.scale_fbr);
                case 'real'
                    [~,h] = contour ( ...
                        space.fbr{1}, ...
                        space.fbr{2}, ...
                        real(fbr) / obj.scale_fbr);
                case 'imag'
                    [~,h] = contour ( ...
                        space.fbr{1}, ...
                        space.fbr{2}, ...
                        imag(fbr) / obj.scale_fbr);
            end
            
        end
        set(h, ...
            'LevelList', linspace(cmin, cmax,  obj.cnt_nlev(1) ), ...
            'LineStyle', obj.patterns{1}, ...
            'LineWidth', obj.l_thin, ...
            'LineColor', obj.colors(m,:))
    end
    
end

if ~obj.range
    axis ( [ space.dof{1}.fbr_min space.dof{1}.fbr_max ...
            space.dof{2}.fbr_min space.dof{2}.fbr_max ] )
else
    axis( [ obj.x_min obj.x_max obj.y_min obj.y_max ] )
end

 
set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)


title ( {info.header1;info.header2} )
if ~obj.marginals
    xlabel ( [ 'P_{', space.dof{1}.label, '}' ] )
end
ylabel ( [ 'P_{', space.dof{2}.label, '}' ] )

