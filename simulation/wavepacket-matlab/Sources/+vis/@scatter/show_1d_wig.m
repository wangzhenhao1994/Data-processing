%---------------------------------------------------------------
% Scatter plots of phase space densities sampled by trajectories
% optionally with expectation values of position / momentum
% optionally with total energy surface
%---------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
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
        if strcmpi(info.system,'Matlab') % why does this not work with Octave?
            h = errorbar ( ...
                expect.pos{1}.cha{m}(obj.mask.ind{m}), ...
                expect.mom{1}.cha{m}(obj.mask.ind{m}), ...
                expect.mom{1}.unc{m}(obj.mask.ind{m}), ...        % lower
                expect.mom{1}.unc{m}(obj.mask.ind{m}), ...        % upper
                expect.pos{1}.unc{m}(obj.mask.ind{m}), ...        % left
                expect.pos{1}.unc{m}(obj.mask.ind{m}) );          % right
            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     obj.colors(m,:), ...
                'LineWidth', obj.l_thick)
        end
    end
    
    % Plot dotted contours of total enery functions
    if obj.energy
        [~,h] = contour ( space.dvr{1} , ...
            -space.fbr{1}/2, ...
            hamilt.tef.wig{m}(:,:) );
        set(h, ...
            'LevelList', linspace( obj.tef_min, obj.tef_max, obj.cnt_nlev ), ...
            'LineStyle', obj.patterns{3}, ...
            'LineColor', obj.colors(m,:), ...
            'LineWidth', obj.l_thin)
    end
    
    % If population not too small, make scatter plot of trajectories
    if expect.pop.cha{m}(step)>expect.min_pop
        
        h = scatter ( state.pos{1}(state.cha==m), state.mom{1}(state.cha==m) );
        set(h, ...
            'Marker',          'o', ...
            'SizeData',        obj.m_extra, ...
            'MarkerFaceColor', obj.colors(m,:), ...
            'MarkerEdgeColor', 'none')
        
        % Surface hopping only
        if isa(state,'sht.generic')
            for  n=1:hamilt.coupling.n_eqs
                
                % Allowed jumps: m <== n: (empty circles)
                allowed = state.Allowed{m,n};
                if ~isempty ( allowed )
                    h = scatter ( state.pos{1}(allowed), state.mom{1}(allowed) );
                    set(h, ...
                        'Marker',          'o', ...
                        'SizeData',        obj.m_extra, ...
                        'MarkerFaceColor', obj.colors(m,:), ...
                        'MarkerEdgeColor', obj.colors(n,:), ...
                        'LineWidth',       obj.l_thick)
                end
                
                % Forbidden jumps: n <== m (crosses)
                forbidden = state.Forbidden{n,m};
                if ~isempty ( forbidden )
                    h = scatter ( state.pos{1}(forbidden), state.mom{1}(forbidden) );
                    set(h, ...
                        'Marker',          'x', ...
                        'SizeData',        obj.m_huger, ...
                        'MarkerEdgeColor', obj.colors(n,:), ...
                        'LineWidth',       obj.l_thick)
                end

                            
            end            
        end
            
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
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
            
            if hamilt.nip{m}.max(1) < space.dof{1}.dvr_max
                h = line ( [ hamilt.nip{m}.max(1) hamilt.nip{m}.max(1)], ...
                    [space.dof{1}.fbr_min/2 space.dof{1}.fbr_max/2] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
        end
    end
end

% Quiver plot of momentum jumps (only for adiabatic two channel surface hopping)
if isa(state,'sht.generic')  ...
        && hamilt.coupling.n_eqs==2  ...
        && strcmpi(hamilt.coupling.represent,'adi')
    
    % Copied from +plot/wigner.m
    kinetic = space.dof{1}.kin;
    kinetic = kinetic(space.dof{1}.n_pts/4+1 : 3*space.dof{1}.n_pts/4);
    kinetic = [kinetic kinetic]';
    kinetic = reshape ( kinetic, space.dof{1}.n_pts, 1 );
    pot_1 = repmat ( hamilt.pot{1,1}.dvr, 1, space.dof{1}.n_pts)';
    pot_2 = repmat ( hamilt.pot{2,2}.dvr, 1, space.dof{1}.n_pts)';

    % Conservation of energy when jumping upwards 2 <== 1
    kin_1 = repmat (        kinetic     , 1, space.dof{1}.n_pts);
    kin_2 = kin_1 - (pot_2 - pot_1);
    mom_1 = sqrt (2*kin_1*space.dof{1}.mass);
    mom_2 = sqrt (2*kin_2*space.dof{1}.mass); 
    
    % Quiver plot of momentum jumps
    if obj.quiver
        delta = mom_2 - mom_1;
        delta(imag(delta)~=0) = 0;
        delta(end/2+1:end,:) = - delta(end/2+1:end,:);
        h = quiver ( space.dvr{1}(1:4:end) , ...
            -space.fbr{1}(1:4:end)/2, ...
            zeros(space.dof{1}.n_pts/4), ...
            delta(1:4:end,1:4:end) );
        set(h, ...
            'Color', obj.colors(1,:), ...
            'AutoScale', 'on')
    end
    
    % Contour separating allowed from forbidden region
    [~,h] = contour ( space.dvr{1} , ...
        -space.fbr{1}/2, ...
        kin_2 );
    set(h, ...
        'LineColor', [1 1 1]/2, ...
        'LevelList', 0, ...
        'LineWidth', obj.l_thick)
    
    % Conservation of energy when jumping downwards: 1 <== 2
    if obj.quiver
        kin_2 = repmat (        kinetic     , 1, space.dof{1}.n_pts);
        kin_1 = kin_2 + (pot_2 - pot_1);
        mom_1 = sqrt (2*kin_1*space.dof{1}.mass);
        mom_2 = sqrt (2*kin_2*space.dof{1}.mass);
        delta = mom_1 - mom_2;
        delta(end/2+1:end,:) = - delta(end/2+1:end,:);
        h = quiver ( space.dvr{1}(1:4:end) , ...
            -space.fbr{1}(1:4:end)/2, ...
            zeros(space.dof{1}.n_pts/4), ...
            delta(1:4:end,1:4:end) );
        set(h, ...
            'Color', obj.colors(2,:), ...
            'AutoScale', 'on')
    end

end


end


