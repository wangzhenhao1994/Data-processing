%------------------------------------------------------------
% Plot 2-dim position densities and potential energy surfaces
%------------------------------------------------------------
function show_2d_dvr ( obj, state, step )
global expect hamilt info plots space

% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs
    
    % Plot "mean trajectory" as solid curve in position space
    % Expectation value of R_1 vs. R_2
    if obj.expect
        if expect.pop.cha{m}(step)>expect.min_pop
            plots.mask.ind{m} = find ( expect.pop.cha{m} > expect.min_pop );
            h = plot ( expect.pos{1}.cha{m}(plots.mask.ind{m}), ...
                expect.pos{2}.cha{m}(plots.mask.ind{m}) );

            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'LineWidth', obj.l_thick, ...
                'Color',     obj.colors(m,:))

        end
    end
       
    % Plot dotted contours of potential energy functions
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
        
        if isa(state,'traj') % this includes sub-classes for SHT
            
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
                
        elseif isa(state,'wave')
            
            switch lower(obj.complex)
                
                case 'abs2'
                    [~,h] = contour ( ...
                        space.dvr{1}, ...
                        space.dvr{2}, ...
                        abs(state.dvr{m}).^2 / obj.scale_dvr);
                case 'real'
                    [~,h] = contour ( ...
                        space.dvr{1}, ...
                        space.dvr{2}, ...
                        real(state.dvr{m}) / obj.scale_dvr);
                case 'imag'
                    [~,h] = contour ( ...
                        space.dvr{1}, ...
                        space.dvr{2}, ...
                        imag(state.dvr{m}) / obj.scale_dvr);
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
if ~obj.marginals
    xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
else
    if strcmpi(info.system,'Matlab')
        set(gca, 'XTickLabels', {})
    elseif strcmpi(info.system,'Octave')
        set(gca, 'XTickLabel', {})
    end
end
ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )

% Absorbing boundary conditions
if isa(state,'wave') && isfield (hamilt,'nip')
    for  m=1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.nip{m}.dvr)
            if hamilt.nip{m}.min(1) > space.dof{1}.dvr_min % Left
                h = line ( [ hamilt.nip{m}.min(1) hamilt.nip{m}.min(1) ], ...
                    [ space.dof{2}.dvr_min space.dof{2}.dvr_max ] );

                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'LineWidth', obj.l_thin, ...
                    'Color',     obj.colors(m,:))

            end
            
            if hamilt.nip{m}.max(1) < space.dof{1}.dvr_max % Right
                h = line ( [ hamilt.nip{m}.max(1) hamilt.nip{m}.max(1) ], ...
                    [ space.dof{2}.dvr_min space.dof{2}.dvr_max ] );

                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'LineWidth', obj.l_thin, ...
                    'Color',     obj.colors(m,:))

            end
            
            if hamilt.nip{m}.min(2) > space.dof{2}.dvr_min % Lower
                h = line ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ], ...
                    [ hamilt.nip{m}.min(2) hamilt.nip{m}.min(2) ] );
                
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'LineWidth', obj.l_thin, ...
                    'Color',     obj.colors(m,:))

            end
            
            if hamilt.nip{m}.max(2) < space.dof{2}.dvr_max % Upper
                h = line ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ], ...
                    [ hamilt.nip{m}.max(2) hamilt.nip{m}.max(2) ] );

                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'LineWidth', obj.l_thin, ...
                    'Color',     obj.colors(m,:))

            end
        end
    end
end

