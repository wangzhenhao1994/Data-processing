%--------------------------------------------------------------------------
%
% Comparison of 'balanced truncation' and/or 'singular perturbation theory'
% and/or 'H2 optimal model reduction' versus full dimensionality dynamics.
% For equations of motions as determined  by input variable "eom" which can
% be 'tdse'|'lvne' in which case we deal with 'ket'|'rho' objects, resp. .
% The input vectors "truncate", "singular" and/or "H2reduce" 
% contain a variable number of approximation degrees to be
% simulated; if one of those vectors is left empty on input, the 
% respective simulation method will be skipped. Graphical output will
% be displayed in figure specified by integer "fig".
% 
% switchoff=0: Do complete simulations
% switchoff>0: Skip calculation of bound states and matrix elements 
% switchoff>1: Skip calculation of A, B, N, C system matrices
% switchoff>2: Skip balancing transformation (BT method only)
%--------------------------------------------------------------------------


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2013-16 Burkhard Schmidt, FU Berlin

function qm_BTversusH2 (eom, fig, truncate, singular, H2reduce, switchoff)

switch eom
    case 'tdse'
        choice = 'ket';
    case 'lvne'
        choice = 'rho';
    otherwise
        prt.error ('Wrong choice for equation of motion')
end

%% Initializes general information and sets up log files.
prt.init (mfilename('fullpath'));

% Initialize: Bound states, matrix elements, full dimensional, balancing
if switchoff<1 % Not necessary if energy/dipole matrices from external sources
    qm_setup('wave'); qm_init; qm_bound;  qm_cleanup; 
    qm_setup('wave'); qm_init; qm_matrix; qm_cleanup; 
end

if switchoff<2 % Not necessary if ABNCD from external sources
    qm_setup(choice); qm_init; qm_abncd; qm_cleanup;
end

%% Control problem in full dimensionality as reference solution
qm_setup(choice); qm_init; qm_propa; qm_cleanup;
global info state time
reference = state.y;

ndim = length(state.x);
for j =1:size(state.y,2)
    mylegend {j,1} = strcat(int2str(ndim), ' (', state.y_label{j}, ')');
end

info.program = 'qm_BTversusH2';

% Open new figure
figure (fig); clf
thisplot = vis.styles; % Construct object
show_logo (thisplot)

%% Simple truncation
if ~isempty (truncate)
    subplot(1,3,1)
    hold on
    set ( gca, ...
        'LineWidth',     thisplot.l_thick, ...
        'FontName',      thisplot.f_name,  ...
        'FontSize',      thisplot.f_large, ...
        'FontWeight',    thisplot.f_heavy )
    for len=1:size(state.y,2)
        plot (time.steps.m_grid,reference(:,len)-reference(1,len),...
            'LineStyle', thisplot.patterns{1}, ...
            'LineWidth', thisplot.l_thin,...
            'Color',thisplot.colors(len,:))
    end
    
    % Balancing transformation
    if switchoff<3
        qm_setup (choice); qm_init; qm_balance;
        switchoff = 3;
    end
    
    %  Loop: different truncations
    for k=1:length(truncate)
        
        % Perform truncation and run in reduced dimensionality
        qm_setup (choice); qm_init; qm_truncate('t',truncate(k));
        qm_setup (choice); qm_init(['t' int2str(truncate(k))]); qm_propa; 
        qm_cleanup;
        
        % Plot resulting populations
        figure(fig);
        subplot(1,3,1)
        global reduce state time
        for len=1:size(state.y,2)
            plot(time.steps.m_grid,real(state.y(:,len)-state.y(1,len)),...
                'LineStyle', thisplot.patterns{k+1}, ...
                'LineWidth', thisplot.l_thin,...
                'Color',thisplot.colors(len,:))
        end
        
        for j =1:length (state.y_label)
            mylegend {size(state.y,2)*k+j,1} = strcat(int2str(truncate(k)), ' (', state.y_label{j}, ')');
        end
        
    end
    hold off
    legend (mylegend, 'Location','NorthEast', 'FontSize',thisplot.f_small)
    xlabel ('time t')
    ylabel ('observables y(t)-y(0)')
    switch lower(reduce.balance.A_stable)
        case 'ssu'
            title  ({'Balancing & Simple truncatation',['SSU stabilization: ' int2str(reduce.balance.A_split) ', B/N scaling: ' num2str(reduce.balance.BN_scale)],['ACF couple: ',int2str(reduce.balance.acf_couple),', Transform: ',upper(reduce.balance.transform)]})
        case 'evs'
            title  ({'Balancing & Simple truncatation',['EVS stabilization: ' num2str(reduce.balance.A_shift) ', B/N scaling: ' num2str(reduce.balance.BN_scale)],['Transform: ',upper(reduce.balance.transform)]})
    end
    drawnow
    
end

%% Singular perturbation
if ~isempty (singular)
    subplot(1,3,2)
    hold on
    set ( gca, ...
        'LineWidth',     thisplot.l_thick, ...
        'FontName',      thisplot.f_name,  ...
        'FontSize',      thisplot.f_large, ...
        'FontWeight',    thisplot.f_heavy )
    for len=1:size(state.y,2)
        plot (time.steps.m_grid,reference(:,len)-reference(1,len),...
            'LineStyle', thisplot.patterns{1}, ...
            'LineWidth', thisplot.l_thin,...
            'Color',thisplot.colors(len,:))
    end
    
    % Balancing transformation
    if switchoff<3
        qm_setup (choice); qm_init; qm_balance;
    end
    
    %  Loop: different truncations
    for k=1:length(singular)
        
        % Perform singular perturbation and run in reduced dimensionality
        qm_setup (choice); qm_init; qm_truncate('s',singular(k));
        qm_setup (choice); qm_init(['s' int2str(singular(k))]); qm_propa;
        qm_cleanup;
        
        % Plot resulting populations
        figure(fig);
        subplot(1,3,2)
        global reduce state time
        for len=1:size(state.y,2)
            plot(time.steps.m_grid,real(state.y(:,len)-state.y(1,len)),...
                'LineStyle', thisplot.patterns{k+1}, ...
                'LineWidth', thisplot.l_thin,...
                'Color',thisplot.colors(len,:))
        end
        
        for j =1:length (state.y_label)
            mylegend {size(state.y,2)*k+j,1} = strcat(int2str(singular(k)), ' (', state.y_label{j}, ')');
        end
        
    end
    hold off
    legend (mylegend, 'Location','NorthEast', 'FontSize',thisplot.f_small)
    xlabel ('time t')
    ylabel ('observables y(t)-y(0)')
    switch lower(reduce.balance.A_stable)
        case 'ssu'
            title  ({'Balancing & Singular perturbation',['SSU stabilization: ' int2str(reduce.balance.A_split) ', B/N scaling: ' num2str(reduce.balance.BN_scale)],['ACF couple: ',int2str(reduce.balance.acf_couple),', Transform: ',upper(reduce.balance.transform)]})
        case 'evs'
            title  ({'Balancing & Singular perturbation',['EVS stabilization: ' num2str(reduce.balance.A_shift) ', B/N scaling: ' num2str(reduce.balance.BN_scale)],['Transform: ',upper(reduce.balance.transform)]})
    end
    drawnow
    
end

%% H2 optimal model reduction
if ~isempty (H2reduce)
    figure(fig)
    subplot(1,3,3)
    hold on
    set ( gca, ...
        'LineWidth',     thisplot.l_thick, ...
        'FontName',      thisplot.f_name,  ...
        'FontSize',      thisplot.f_large, ...
        'FontWeight',    thisplot.f_heavy )
    
    for len=1:size(state.y,2)
        plot (time.steps.m_grid,real(reference(:,len)-reference(1,len)),...
            'LineStyle', thisplot.patterns{1}, ...
            'LineWidth', thisplot.l_thin,...
            'Color',thisplot.colors(len,:))
    end
    
    for k=1:length(H2reduce)
        
        % Perform dimension reduction and run in reduced dimensionality
        qm_setup(choice); qm_init; qm_H2model (H2reduce(k));
        qm_setup(choice); qm_init(['h',int2str(H2reduce(k))]); qm_propa;
        qm_cleanup;
        
        % Plot resulting populations
        figure(fig);
        subplot(1,3,3)
        global reduce state time
        for len=1:size(state.y,2)
            plot(time.steps.m_grid, real(state.y(:,len)-state.y(1,len)),...
                'LineStyle', thisplot.patterns{k+1}, ...
                'LineWidth', thisplot.l_thin,...
                'Color',thisplot.colors(len,:))
        end
        
        for j =1:length (state.y_label)
            mylegend {size(state.y,2)*k+j,1} = strcat(int2str(H2reduce(k)), ' (', state.y_label{j}, ')');
        end
        
    end
    hold off
    legend (mylegend, 'Location','NorthEast', 'FontSize',thisplot.f_small)
    xlabel ('time t')
    ylabel ('observables y(t)-y(0)')
    switch lower(reduce.H2model.A_stable)
        case 'ssu'
            title  ({'H2 model reduction',['SSU stabilization: ' int2str(reduce.H2model.A_split) ', B/N scaling: ' num2str(reduce.H2model.BN_scale)],['Tolerance: ',num2str(reduce.H2model.conv_tol),', Max. iterations: ',int2str(reduce.H2model.max_iter)]})
        case 'evs'
            title  ({'H2 model reduction',['EVS stabilization: ' num2str(reduce.H2model.A_shift) ', B/N scaling: ' num2str(reduce.H2model.BN_scale)],['Tolerance: ',num2str(reduce.H2model.conv_tol),', Max. iterations: ',int2str(reduce.H2model.max_iter)]})
    end
    
end

% Output clock/date/time
prt.clock;

% Save figure
saveas (gca, int2str(fig), 'fig')

end


