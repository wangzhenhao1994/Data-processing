%--------------------------------------------------------------------------
% Lower part: flux density
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2013-2017 Burkhard Schmidt
%
% see the README file for license details.

function show_1d_lower ( obj, state, step )
global expect hamilt space

% Loop over all channels where population not too small
for  m=1:hamilt.coupling.n_eqs
    if ( expect.pop.cha{m}(step) > expect.min_pop )
        
        if isa (state,'traj')
            
            % Stem plot of flux (momenta)
            h = stem ( state.pos{1}(state.cha==m), state.mom{1}(state.cha==m) );
            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     obj.colors(m,:), ...
                'LineWidth', obj.l_thick)
            
        elseif isa (state,'wave')
            
            % Calculate flux density
            psi_mom = momentum(space.dof{1}, state.dvr{m});
            j = real ( conj(state.dvr{m}) .* psi_mom );
            
            % Curve plot of flux density
            h = plot ( space.dvr{1}, j);
            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     obj.colors(m,:), ...
                'LineWidth', obj.l_thick)
            
        end
    end
end

%% Axes, labels, etc
if ~obj.range
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max space.dof{1}.fbr_min/2 space.dof{1}.fbr_max/2 ] )
else
    axis ( [obj.x_min obj.x_max obj.y_min obj.y_max ] )
end

set(gca, ...
    'XAxisLocation', 'bottom', ...
    'LineWidth',     obj.l_thick, ...
    'FontName',      obj.f_name, ...
    'FontSize',      obj.f_large, ...
    'FontWeight',    obj.f_heavy)

xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
ylabel('j(R)' )

end

