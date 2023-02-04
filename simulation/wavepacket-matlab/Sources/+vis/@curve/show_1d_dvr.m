%-------------------------------------------------
% Plot position density and potential energy curve
%-------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-20xx Burkhard Schmidt's group
%               2009 Ulf Lorenz
%
% see the README file for license details.

function show_1d_dvr ( obj, state, step )
global expect hamilt info space


%% Plot potential energy curve (optional)
if obj.energy
    
    e_label = cell(1,hamilt.coupling.n_eqs);
    sub_set = [];
    
    for  m=1:hamilt.coupling.n_eqs
        if obj.energy
            if ~isempty(hamilt.pot{m,m}.dvr)
                h = plot ( space.dvr{1}, hamilt.pot{m,m}.dvr);
            else
                h = plot ( space.dvr{1}, zeros(size(space.dvr{1})));
            end
            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     obj.colors(m,:), ...
                'LineWidth', obj.l_thin)
            
            % Labelling coupled states
            sub_set = [sub_set h];
            if strcmpi (hamilt.coupling.represent,'dia')
                e_label{1,m} = hamilt.coupling.labels{m};
            elseif strcmpi (hamilt.coupling.represent,'adi')
                e_label{1,m} = ['adi-' int2str(m)];
            end
        end
        
    end
    
    legend (sub_set, e_label, 'AutoUpdate', 'off')
    
end

%% Loop over (coupled) channels
for  m=1:hamilt.coupling.n_eqs
    
    % If population not too small
    if ( expect.pop.cha{m}(step) > expect.min_pop )
        
        % Horizontal offset (= kinetic + potential energy)
        if obj.energy
            offset = expect.kin.cha{m}(step) + expect.pot.cha{m}(step);
            scaling = obj.pot_delta;
            y_min = obj.pot_min-obj.pot_delta/10;
            y_max = obj.pot_max+obj.pot_delta/10;
        else
            offset = 0;
            scaling = 1;
            y_min = -0.1;
            y_max = +1.1;
        end
        
        if isa(state,'traj')  % this includes sub-classes for SHT
            x  = space.dof{1}.x_grid;
            n  = space.dof{1}.n_pts;
            dx = space.dof{1}.x_dlt;
            
            if strcmpi(info.system,'Matlab') % histcounts not available in Octave
                [rho,~] = histcounts (state.pos{1}(state.cha==m), x);
                h = plot ( x(1:n-1)+dx/2, rho*scaling/obj.scale_dvr+offset);
            elseif strcmpi(info.system,'Octave') % use of histc is discouraged in Matlab
                rho = histc (state.pos{1}(state.cha==m), x)';
                h = plot ( x(1:n), rho*scaling/obj.scale_dvr+offset);
            end
            
           set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     obj.colors(m,:), ...
                'LineWidth', obj.l_thick)
            
        elseif isa (state,'wave')
            
            % Three options of what to plot as curve(s)
            switch lower(obj.complex)
                case 'abs2'
                    % Plot density as height with phase as color: Map [-pi,pi] into [0,1]
                    rho = abs   ( state.dvr{m} ) .^2;
                    phi = angle ( state.dvr{m} ) / (2*(pi+0.001)) + 1/2;
                    if strcmpi(info.system,'Matlab')
                        vis.styles.plot_color ( space.dvr{1}, ...
                            rho*scaling/obj.scale_dvr, ...
                            phi, ...
                            obj.colors(m,:), ...
                            obj.l_extra, ...
                            offset, ...
                            0 )
                    elseif  strcmpi(info.system,'Octave')
                        h = plot ( space.dvr{1}, rho*scaling/obj.scale_dvr + offset);
                        set(h, ...
                            'LineStyle', obj.patterns{1}, ...
                            'LineWidth', obj.l_thick, ...
                            'Color',     obj.colors(m,:) )
                    end
                    
                case 'real'
                    % Plot real part of wavefunction
                    wav = real ( state.dvr{m} );
                    h = plot ( space.dvr{1}, ...
                        wav*scaling/sqrt(obj.scale_dvr)+offset );
                    set(h, ...
                        'LineStyle', obj.patterns{1}, ...
                        'Color',     obj.colors(m,:), ...
                        'LineWidth', obj.l_thick)

                case 'imag'
                    % Plot real part of wavefunction
                    wav = imag ( state.dvr{m} );
                    h = plot ( space.dvr{1}, ...
                        wav*scaling/sqrt(obj.scale_dvr)+offset );
                    set(h, ...
                        'LineStyle', obj.patterns{1}, ...
                        'Color',     obj.colors(m,:), ...
                        'LineWidth', obj.l_thick)

                otherwise
                    prt.error ('Wrong choice of what to plot')
            end
            
        end
        
        % Plot horizontal base line
        if obj.energy
            h = line ( [space.dof{1}.dvr_min space.dof{1}.dvr_max],[1 1]*offset );

            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     obj.colors(m,:), ...
                'LineWidth', obj.l_thin)

        end
        
    end
    
end

%% Axes, labels, etc
if ~obj.range
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max y_min y_max ] )
else
    if ~isempty(obj.y_min) && ~isempty(obj.y_max)
        axis ( [ obj.x_min obj.x_max obj.y_min obj.y_max ] )
    else
        axis ( [ obj.x_min obj.x_max y_min y_max ] )
    end
end

set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)

title ( {info.header1;info.header2} )
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )

if obj.energy
    
    if hamilt.coupling.n_eqs==1
        ylabel ( 'V(R)' )
    else
        if strcmpi ( hamilt.coupling.represent,'adi')
            ylabel ( 'V_{adi}(R)' )
        elseif strcmpi ( hamilt.coupling.represent,'dia')
            ylabel ( 'V_{dia}(R)' )
        end
    end
    
else
    
    if hamilt.coupling.n_eqs==1
        ylabel ( '\rho(R)' )
    else
        if strcmpi ( hamilt.coupling.represent,'adi')
            ylabel ( '\rho_{adi}(R)' )
        elseif strcmpi ( hamilt.coupling.represent,'dia')
            ylabel ( '\rho_{dia}(R)' )
        end
    end
    
end

% Negative imaginary potential (as absorbing boundary condition)
if isa (state,'wave') && isfield (hamilt,'nip')
    for  m=1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.nip{m}.dvr)
            if hamilt.nip{m}.min(1) > space.dof{1}.dvr_min
                h = line ( [hamilt.nip{m}.min(1) hamilt.nip{m}.min(1)], ...
                    [obj.pot_min-obj.pot_delta/10 obj.pot_max+obj.pot_delta/10] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)

            end
            
            if hamilt.nip{m}.max(1) < space.dof{1}.dvr_max
                h = line ( [hamilt.nip{m}.max(1) hamilt.nip{m}.max(1)], ...
                    [obj.pot_min-obj.pot_delta/10 obj.pot_max+obj.pot_delta/10] );
                set(h, ...
                    'LineStyle', obj.patterns{2}, ...
                    'Color',     obj.colors(m,:), ...
                    'LineWidth', obj.l_thin)
            end
        end
    end
end

end


