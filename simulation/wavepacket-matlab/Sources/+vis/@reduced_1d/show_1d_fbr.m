%--------------------------------------------------------------------------
% Plot reduced density matrices for k-th dimension
% in momentum (P_k) representation (Q/C only)
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20xy Burkhard Schmidt's group
%
% see the README file for license details.

function show_1d_fbr ( obj, state, step, k )
global expect hamilt info space

if isa(state,'wave')
    prt.error ('Momentum representation not available for wavefunctions')
end

% Initialize calculations of densities and their maxima
rho = cell(hamilt.coupling.n_eqs,1);
if step==1
    obj.rho_max(k) = 0;
end

% Loop over (coupled) densities
for  m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step) > expect.min_pop
        
        y  = space.dof{k}.p_grid;
        n  = space.dof{k}.n_pts;
        dy = space.dof{k}.p_dlt;
        if strcmpi(info.system,'Matlab') % histcounts not available in Octave
            [rho{m},~] = histcounts ( state.mom{k}(state.cha==m), y);
        elseif strcmpi(info.system,'Octave') % use of histc is discouraged in Matlab
            rho{m} = histc ( state.mom{k}(state.cha==m), y);
        end
        
        % Find maximal density at first step
        if step==1
            obj.rho_max(k) = max ( obj.rho_max(k), max(abs(rho{m})) );
        end
        
    end
end

% Loop over (coupled) densities
for  m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step) > expect.min_pop

        % Create curve plots
        if strcmpi(info.system,'Matlab')
            h = plot ( y(1:n-1)+dy/2, rho{m}/obj.rho_max(k) );
        elseif strcmpi(info.system,'Octave')
            h = plot ( y(1:n), rho{m}/obj.rho_max(k) );
        end
        set(h, ...
            'LineStyle', obj.patterns{1}, ...
            'Color',     obj.colors(m,:), ...
            'LineWidth', obj.l_thick)
    end
end

% Optionally setting plot ranges "manually"
if ~obj.range
    axis ( [ space.dof{k}.fbr_min space.dof{k}.fbr_max -0.1 +1.1 ] )
else
    axis ( [ obj.x_min            obj.x_max            -0.1 +1.1 ] )
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
xlabel ( [ 'P_{', space.dof{k}.label, '}' ] )
if hamilt.coupling.n_eqs==1
    ylabel ( '\rho(P)' )
else
    if strcmpi ( hamilt.coupling.represent,'adi' )
        ylabel ( '\rho_{adi}(P)' )
    elseif strcmpi ( hamilt.coupling.represent,'dia' )
        ylabel ( '\rho_{dia}(P)' )
    end
end

end


