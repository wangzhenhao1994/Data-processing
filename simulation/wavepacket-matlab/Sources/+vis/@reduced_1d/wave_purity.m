
%--------------------------------------------------------------------
% Plot purity measure for k-th dimension vs. time
%--------------------------------------------------------------------
function wave_purity ( obj, psi, step, k )
global expect info time space hamilt

% Calculate tr(rho^2)
expect.pur{k}.tot(step) = 0;
for m = 1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step) > expect.min_pop
        expect.pur{k}.cha{m}(step) = abs(trace(psi.redu{m,k}^2));
        expect.pur{k}.tot   (step) = expect.pur{k}.tot(step) + expect.pur{k}.cha{m}(step) / expect.pop.cha{m}(step)^2;
    end
end

% Plot evolving curve
mask = find ( expect.pur{k}.tot > eps );
h = plot ( time.steps.m_grid (mask), ...
    expect.pur{k}.tot(mask) );
set(h, ...
    'LineStyle', obj.patterns{1}, ...
    'Color',     'black', ...
    'LineWidth', obj.l_thick)

% Line styles and fonts
set(gca, ...
	'LineWidth',  obj.l_thick, ...
	'FontName',   obj.f_name, ...
	'FontSize',   obj.f_large, ...
	'FontWeight', obj.f_heavy)

% Axes range and labels
axis ( [ 0 time.steps.t_total -0.1 1.1 ] )
if strcmpi(info.program,'qm_propa')
    xlabel ('t')
elseif strcmpi(info.program,'qm_bound')
    xlabel ('n')
end
    
switch space.n_dim
    case 1
        ylabel ( 'tr(\rho^2)' )
    case 2
        ylabel ( 'tr(\rho_1^2) = tr(\rho_2^2)' )
    otherwise
        ylabel ( ['tr(\rho_' int2str(k) '^2)'] )
end   

