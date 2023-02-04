%------------------------------------------------------------
% Plot reduced density matrices for k-th dimension
% in position (R_k) representation
%------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20xy Burkhard Schmidt's group
%
% see the README file for license details.

function show_1d_dvr ( obj, state, step, k )       
global info space expect hamilt

% Initialize calculations of densities and their maxima
rho = cell(hamilt.coupling.n_eqs,1);
if step==1
    obj.rho_max(k) = 0;
end

% Main switch: Q/M versus Q/C
if isa(state,'traj')  % this includes sub-classes for SHT
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Get 1-dim histograms from trajectory bundles
            x  = space.dof{k}.x_grid;
            n  = space.dof{k}.n_pts;
            dx = space.dof{k}.x_dlt;
            if strcmpi(info.system,'Matlab') % histcounts not available in Octave
                [rho{m},~] = histcounts ( state.pos{k}(state.cha==m), x );
            elseif strcmpi(info.system,'Octave') % use of histc is discouraged in Matlab
                rho{m} = histc ( state.pos{k}(state.cha==m), x );
            end
            
            % Find maximal density at first step
            if step==1
                obj.rho_max(k) = max ( obj.rho_max(k), max(abs(rho{m})) );
            end
            
        end
    end
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Create curve plots
            if strcmpi(info.system,'Matlab') 
                h = plot ( x(1:n-1)+dx/2, rho{m}/obj.rho_max(k) );
            elseif strcmpi(info.system,'Octave')
                h = plot ( x(1:n), rho{m}/obj.rho_max(k) );
            end
            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     obj.colors(m,:), ...
                'LineWidth', obj.l_thick)
        end
    end
    
elseif isa(state,'wave')
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Get diagonal elements of reduced density matrices
            rho{m} = diag(state.redu{m,k});
            
            % Find maximal density at first step
            if step==1
                obj.rho_max(k) = max ( obj.rho_max(k), max(rho{m}) );
            end
            
        end
    end
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Create curve plots
            h = plot ( space.dof{k}.x_grid, rho{m}/obj.rho_max(k) );
            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     obj.colors(m,:), ...
                'LineWidth', obj.l_thick)
        end
    end
    
end

% Optionally setting plot ranges "manually"
if ~obj.range
    axis ( [ space.dof{k}.dvr_min space.dof{k}.dvr_max -0.1 +1.1 ] )
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
xlabel ( ['R_{', space.dof{k}.label, '}'] )
if hamilt.coupling.n_eqs==1
    ylabel ( '\rho(R)' )
else
    if strcmpi ( hamilt.coupling.represent,'adi' )
        ylabel ( '\rho_{adi}(R)' )
    elseif strcmpi ( hamilt.coupling.represent,'dia' )
        ylabel ( '\rho_{dia}(R)' )
    end
end

% Negative imaginary potential (as absorbing boundary conditions)
if isa(state,'wave') && isfield (hamilt,'nip')
    for  m=1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.nip{m}.dvr)
            if hamilt.nip{m}.min(k) > space.dof{k}.dvr_min % Left border
                h = line ( [ hamilt.nip{m}.min(k) hamilt.nip{m}.min(k) ], ...
                    [ -0.1 +1.1 ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
            
            if hamilt.nip{m}.max(k) < space.dof{k}.dvr_max % Right border
                h = line ( [ hamilt.nip{m}.max(k) hamilt.nip{m}.max(k) ], ...
                    [ -0.1 +1.1 ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
        end
    end
end
