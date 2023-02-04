%--------------------------------------------------------------------
% Plot autocorrelation vs. time
%--------------------------------------------------------------------
function correlation (obj)
global info time

% Not too many time steps: plot colored curve
if time.steps.s_number * time.steps.m_number<=1000
    
    % Get modulus of autocorrelation function
    rho = abs   ( time.steps.acf ) .^2;
    
    % Get phase of autocorrelation function: Map interval [-pi,pi] into [0,1]
    phi = angle ( time.steps.acf ) / (2*(pi+0.001)) + 1/2;
    
    % Plot autocorrelation function (colored curve)
    if strcmpi(info.system,'Matlab')
        vis.styles.plot_color ( time.steps.s_grid, rho, phi, [1 1 0]/2, obj.l_thick, 0, 0 )
    elseif  strcmpi(info.system,'Octave')
        h = plot ( time.steps.s_grid, rho );
        set(h, ...
            'LineStyle', obj.patterns{1}, ...
            'Color',     'black', ...
            'LineWidth', obj.l_thick)
        
    end
    
    % Labelling the ordinate
    my_label = 'C(t)';
    
    % Too many time steps: plot black curve showing abs(C) instead
else

    h = plot ( time.steps.s_grid, abs(time.steps.acf) );
    set(h, ...
        'LineStyle', obj.patterns{1}, ...
        'Color',     'black', ...
        'LineWidth', obj.l_thick)
    
    % Labelling the ordinate
    my_label = 'autocorrelation';
end

% Axes and labels
axis ( [ 0 time.steps.t_total -0.1 1.1 ] )
set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)

xlabel ('t')
ylabel (my_label)

end
