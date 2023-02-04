%----------------------------------------------------------
% Contour plots of phase space densities
% optionally with expectation values of position / momentum
% optionally with total energy surface
% optionally with boundaries of negative imaginary potential
%----------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function show_1d_wig ( obj, state, step )
global expect hamilt info space

% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs
    
    % Plot "mean trajectory" as solid curve with error bars in phase space:
    % Expectation value of momentum versus expectation value of position
    if obj.expect
        obj.mask.ind{m} = find ( expect.pop.cha {m}(:) > expect.min_pop );
        h = errorbar ( ...
            expect.pos{1}.cha{m}(obj.mask.ind{m}), ...
            expect.mom{1}.cha{m}(obj.mask.ind{m}), ...
            expect.mom{1}.unc{m}(obj.mask.ind{m}), ... % lower
            expect.mom{1}.unc{m}(obj.mask.ind{m}), ... % upper
            expect.pos{1}.unc{m}(obj.mask.ind{m}), ... % left
            expect.pos{1}.unc{m}(obj.mask.ind{m}) ); % right
        set(h, ...
            'LineStyle', obj.patterns{1}, ...
            'LineWidth', obj.l_thick, ...
            'Color',     obj.colors(m,:))
    end
    
    % Plot dotted contours of total enery functions
    if obj.energy
        [~,h] = contour ( space.dvr{1} , ...
            -space.fbr{1}/2, ...
            hamilt.tef.wig{m}(:,:) );
        
        set(h, ...
            'LevelList', linspace( obj.tef_min, obj.tef_max, obj.cnt_nlev(2) ), ...
            'LineStyle', obj.patterns{3}, ...
            'LineWidth', obj.l_thin, ...
            'LineColor', obj.colors(m,:) )
        
    end
    
    % If population not too small, plot solid contours of Wigner functions
    if expect.pop.cha{m}(step)>expect.min_pop
        
        if isa(state,'traj') % this includes sub-classes for SHT
            if obj.cnt_levels
                cmin = obj.cnt_min;
                cmax = obj.cnt_max;
            else
                cmin = 1 / obj.cnt_nlev(1);
                cmax = 1;
            end
            
            x  = space.dof{1}.x_grid;
            y  = space.dof{1}.p_grid;
            n  = space.dof{1}.n_pts;
            dx = space.dof{1}.x_dlt;
            dy = space.dof{1}.p_dlt;
            
            if strcmpi(info.system,'Matlab') % histcounts2 not available in Octave
                [c,~] = histcounts2 ( state.pos{1}(state.cha==m), state.mom{1}(state.cha==m), x, y );
                [~,h] = contour ( x(1:n-1)+dx/2, y(1:n-1)+dy/2, c'/obj.scale_wig );
            elseif strcmpi(info.system,'Octave') % use of hist3 is discouraged in Matlab
                [c,~] = hist3 ( [state.pos{1}(state.cha==m), state.mom{1}(state.cha==m)], 'Ctrs',{x y} );
                [~,h] = contour ( x(1:n), y(1:n), c'/obj.scale_wig );
            end
            
        elseif isa(state,'wave')
            if obj.cnt_levels
                cmin = obj.cnt_min;
                cmax = obj.cnt_max;
            else
                cmin = -state.wig_max;
                cmax = +state.wig_max;
            end
            
            [~,h] = contour ( space.dvr{1},    ...
                -space.fbr{1}/2, ...
                real(state.wig{m}) );
        end
        
        set(h, ...
            'LevelList', linspace(cmin, cmax, obj.cnt_nlev(1) ), ... % use even number of contours to avoid zero! 
            'LineStyle', obj.patterns{1}, ...
            'LineWidth', obj.l_thin, ...
            'LineColor', obj.colors(m,:) )

    end
    
end

% Optionally setting plot ranges "manually"
if ~obj.range
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
        space.dof{1}.fbr_min/2 space.dof{1}.fbr_max/2 ] )
else
    axis ( [obj.x_min obj.x_max obj.y_min obj.y_max ] )
end

set(gca, ...
    'XAxisLocation', 'top', ...
    'LineWidth',     obj.l_thick, ...
    'FontName',      obj.f_name, ...
    'FontSize',      obj.f_large, ...
    'FontWeight',    obj.f_heavy)

xlabel ( {info.header1;info.header2} )
ylabel ( [ 'P_{', space.dof{1}.label, '}' ] )

% Negative imaginary potential (as absorbing boundary conditions)
if isa(state,'wave') && isfield (hamilt,'nip')
    for  m=1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.nip{m}.dvr)
            if hamilt.nip{m}.min(1) > space.dof{1}.dvr_min
                h = line ( [ hamilt.nip{m}.min(1) hamilt.nip{m}.min(1)], ...
                    [space.dof{1}.fbr_min/2 space.dof{1}.fbr_max/2] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'LineWidth', obj.l_thin, ...
                    'Color',     obj.colors(m,:))
            end
            
            if hamilt.nip{m}.max(1) < space.dof{1}.dvr_max
                h = line ( [ hamilt.nip{m}.max(1) hamilt.nip{m}.max(1)], ...
                    [space.dof{1}.fbr_min/2 space.dof{1}.fbr_max/2] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'LineWidth', obj.l_thin, ...
                    'Color',     obj.colors(m,:))
            end
        end
    end
end

end


