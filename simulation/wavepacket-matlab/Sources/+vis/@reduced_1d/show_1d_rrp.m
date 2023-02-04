%-------------------------------------------------------------
% Contour plots of reduced density matrices for k-th dimension
% in position (R_k,R_k') representation (Q/M only)
%-------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20xy Burkhard Schmidt's group
%
% see the README file for license details.

function show_1d_rrp ( obj, state, step, k )
global info space expect hamilt

% Initialize calculations of maximal densities
if step==1
    obj.rho_max(k)=0;
end

if ~isa(state,'wave')
    prt.error ('R-R'' representation available only for wavefunctions')
end

% Loop over (coupled) densities
for m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step)>expect.min_pop
        
        % Find maximal density at first step
        if step==1
            obj.rho_max(k) = max ( obj.rho_max(k), max(max(abs(state.redu{m,k}))) );
        end
        
    end
end
        
% Loop over (coupled) densities
for m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step)>expect.min_pop

        % Create contour plots
        [~,h] = contour ( ...
            space.dof{k}.x_grid, ...
            space.dof{k}.x_grid, ...
            abs(state.redu{m,k}));
        set(h, ...
            'LevelList', linspace(0, obj.rho_max(k), obj.cnt_nlev), ...   % use even number of contours to avoid zero!
            'LineStyle', obj.patterns{1}, ...
            'LineColor', obj.colors(m,:), ...
            'LineWidth', obj.l_thin)
    end
end

% Dotted line along diagonal
h = line ( [ space.dof{k}.dvr_min space.dof{k}.dvr_max  ], ...
    [ space.dof{k}.dvr_min space.dof{k}.dvr_max ] );
set(h, ...
    'LineStyle', obj.patterns{3}, ...
    'Color',     'black', ...
    'LineWidth', obj.l_thin)

% Optionally setting plot ranges "manually"
if ~obj.range
    axis ( [ space.dof{k}.dvr_min space.dof{k}.dvr_max ...
        space.dof{k}.dvr_min space.dof{k}.dvr_max ] )
else
    axis ( [obj.x_min obj.x_max obj.y_min obj.y_max ] )
end

set(gca, ...
	'LineWidth',  obj.l_thick, ...
	'FontName',   obj.f_name, ...
	'FontSize',   obj.f_large, ...
	'FontWeight', obj.f_heavy)

% Place the header approximately in the middle
if k==2
    title ( {info.header1;info.header2} )
end

xlabel ( ['R_{', space.dof{k}.label, '}'] )
ylabel ( ['R_{', space.dof{k}.label, '}'] )

% Negative imaginary potential (as absorbing boundary conditions)
if isfield (hamilt,'nip')
    for  m=1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.nip{m}.dvr)
            if hamilt.nip{m}.min(k) > space.dof{k}.dvr_min % Left/lower border
                h = line ( [ hamilt.nip{m}.min(k) hamilt.nip{m}.min(k) ], ...
                    [ space.dof{k}.dvr_min space.dof{k}.dvr_max ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)

                h = line ( [ space.dof{k}.dvr_min space.dof{k}.dvr_max ], ...
                    [ hamilt.nip{m}.min(k) hamilt.nip{m}.min(k) ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
            
            if hamilt.nip{m}.max(k) < space.dof{k}.dvr_max % Right/upper border
                h = line ( [ hamilt.nip{m}.max(k) hamilt.nip{m}.max(k) ], ...
                    [ space.dof{k}.dvr_min space.dof{k}.dvr_max ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)

                h = line ( [ space.dof{k}.dvr_min space.dof{k}.dvr_max ], ...
                    [ hamilt.nip{m}.max(k) hamilt.nip{m}.max(k) ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
        end
    end
end
