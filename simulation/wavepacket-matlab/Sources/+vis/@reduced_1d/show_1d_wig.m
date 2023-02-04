%-------------------------------------------------------------
% Contour plots of reduced density matrices for k-th dimension 
% in (Wigner) phase space (R_k, P_k) representation
%-------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20xy Burkhard Schmidt's group
%
% see the README file for license details.

function show_1d_wig ( obj, state, step, k )
global info space expect hamilt

% Initialize calculations of densities and their maxima
rho = cell(hamilt.coupling.n_eqs,1);
if step==1
    obj.rho_max(k)=0;
end

% Plot "mean trajectory" as solid curve in phase space
if obj.expect
    for  m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            plots.mask.ind{m} = find ( expect.pop.cha{m} > expect.min_pop );
            h = plot ( expect.pos{k}.cha{m}(plots.mask.ind{m}), ...
                expect.mom{k}.cha{m}(plots.mask.ind{m}) );
            
            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'LineWidth', obj.l_thick, ...
                'Color',     obj.colors(m,:))
        end
    end
end

% Main switch: Q/M versus Q/C
if isa(state,'traj')  % this includes sub-classes for SHT
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Get 2-dim histograms from trajectory bundles
            x  = space.dof{k}.x_grid;
            y  = space.dof{k}.p_grid;
            n  = space.dof{k}.n_pts;
            dx = space.dof{k}.x_dlt;
            dy = space.dof{k}.p_dlt;
            
            if strcmpi(info.system,'Matlab') % histcounts2 not available in Octave
                [rho{m},~] = histcounts2 ( state.pos{k}(state.cha==m), state.mom{k}(state.cha==m), x, y );
            elseif strcmpi(info.system,'Octave') % use of hist3 is discouraged in Matlab
                [rho{m},~] = hist3 ( [state.pos{k}(state.cha==m), state.mom{k}(state.cha==m)], 'Ctrs',{x y} );
            end
            
            % Find maximal density at first step
            if step==1
                obj.rho_max(k) = max ( obj.rho_max(k), max(max(abs(rho{m}))) );
            end
            
        end
    end
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Create contour plots
            if strcmpi(info.system,'Matlab')
                [~,h] = contour ( x(1:n-1)+dx/2, y(1:n-1)+dy/2, rho{m}'/obj.rho_max(k) );
            elseif strcmpi(info.system,'Octave')
                [~,h] = contour ( x(1:n), y(1:n), rho{m}'/obj.rho_max(k) );
            end
            
            set(h, ...
                'LevelList', linspace(0, 1, obj.cnt_nlev(1)), ...   % use even number of contours to avoid zero!
                'LineStyle', obj.patterns{1}, ...
                'LineColor', obj.colors(m,:), ...
                'LineWidth', obj.l_thin)
        end
    end
    
    
elseif isa (state,'wave')
    
    % For FFT grids only
    if ~isa(space.dof{k},'dof.fft')
        prt.error ('Wigner transform not available for non-FFT grids')
    end
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Get Wigner transforms of reduced density matrices
            n = space.dof{k}.n_pts;
            rho{m} = zeros(n);
            for ii=2:n
                off = ii-(n/2+1);
                rho{m}(ii,1+abs(off):n-abs(off)) = diag(state.redu{m,k},2*off);
            end
            rho{m} = fftshift(ifft(fftshift(rho{m},1)),1)/sqrt(2);
            
            % Find maximal density at first step
            if step==1
                obj.rho_max(k) = max ( obj.rho_max(k), max(max(abs(rho{m}))) );
            end
            
        end
    end
    
    % Loop over (coupled) densities
    for m=1:hamilt.coupling.n_eqs
        if expect.pop.cha{m}(step)>expect.min_pop
            
            % Create contour plots
            [~,h] = contour ( ...
                +space.dof{k}.x_grid, ...
                -space.dof{k}.p_grid/2, ...
                real(rho{m}));
            set(h, ...
                'LevelList', linspace(-obj.rho_max(k), obj.rho_max(k), obj.cnt_nlev), ...  % use even number of contours to avoid zero!
                'LineStyle', obj.patterns{1}, ...
                'LineColor', obj.colors(m,:), ...
                'LineWidth', obj.l_thin)
        end
    end
end

% Optionally setting plot ranges "manually"
if ~obj.range
    axis ( [ space.dof{k}.dvr_min space.dof{k}.dvr_max ...
        space.dof{k}.fbr_min/2 space.dof{k}.fbr_max/2 ] )
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
ylabel ( ['P_{', space.dof{k}.label, '}'] )

% Negative imaginary potential (as absorbing boundary conditions)
if isa(state,'wave') && isfield (hamilt,'nip')
    for  m=1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.nip{m}.dvr)
            if hamilt.nip{m}.min(k) > space.dof{k}.dvr_min % Left border
                h = line ( [ hamilt.nip{m}.min(k) hamilt.nip{m}.min(k) ], ...
                    [ space.dof{k}.dvr_min space.dof{k}.dvr_max ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
            
            if hamilt.nip{m}.max(k) < space.dof{k}.dvr_max % Right border
                h = line ( [ hamilt.nip{m}.max(k) hamilt.nip{m}.max(k) ], ...
                    [ space.dof{k}.dvr_min space.dof{k}.dvr_max ] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
        end
    end
end
