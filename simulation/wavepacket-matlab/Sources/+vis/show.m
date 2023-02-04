%--------------------------------------------------------------------------
% Create animated graphics of time evolving densities
% Create animated graphics of time evolving expectation values
% Create graphics of absorption spectrum (FFT of autocorrelation)
%
% If "state" is an object of class "wave": densities from wavefunctions
%
% If "state" is an object of class "traj": densities from trajectories
%                            Similarly for SHT variants: mssh, fssh, sssh
%
% If "state" is an object of class "ket": populations from ket vectors
%
% If "state" is an object of class "rho": populations/coherences from
%                                         density matrices
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2018 Burkhard Schmidt's group
%               2007-2010 Ulf Lorenz
%
% see the README file for license details.

function show ( state, step )
global expect hamilt info plots space time
persistent writerObj

% Fake time axis for the case of only 1 step
if time.steps.t_total==0
    time.steps.t_total=1;
end

% Various initializations of density plots
if step==1 && isfield(plots,'density')
    
    %% Q/C only: Find maxima of binned densities
    if isa(state,'traj')
        if ~isa(plots.density,'vis.reduced_1d') && ~isa(plots.density,'vis.reduced_2d')
            
            dvr_max = 0;
            fbr_max = 0;
            wig_max = 0;
            
            for m=1:hamilt.coupling.n_eqs
                if expect.pop.cha{m}(step)>expect.min_pop
                    
                    switch space.n_dim
                        case 1
                            if strcmpi(info.system,'Matlab') % histcounts[2] not available in Octave
                                [c_dvr,~] = histcounts  (state.pos{1},               space.dof{1}.x_grid                     );
                                [c_fbr,~] = histcounts  (              state.mom{1},                      space.dof{1}.p_grid);
                                [c_wig,~] = histcounts2 (state.pos{1}, state.mom{1}, space.dof{1}.x_grid, space.dof{1}.p_grid);
                            elseif strcmpi(info.system,'Octave') % use of histc/hist3 is discouraged in Matlab
                                c_dvr = histc (state.pos{1},                         space.dof{1}.x_grid                     );
                                c_fbr = histc (               state.mom{1},                              space.dof{1}.p_grid );
                                c_wig = hist3 ([state.pos{1}, state.mom{1}], 'Ctrs',{space.dof{1}.x_grid space.dof{1}.p_grid});
                            end
                        case 2
                            if strcmpi(info.system,'Matlab') % histcounts2 not available in Octave
                                [c_dvr,~] = histcounts2 (state.pos{1}, state.pos{2}, space.dof{1}.x_grid, space.dof{2}.x_grid);
                                [c_fbr,~] = histcounts2 (state.mom{1}, state.mom{2}, space.dof{1}.p_grid, space.dof{2}.p_grid);
                            elseif strcmpi(info.system,'Octave') % use of hist3 is discouraged in Matlab                            
                                [c_dvr,~] = hist3 ([state.pos{1}, state.pos{2}], 'Ctrs',{space.dof{1}.x_grid space.dof{2}.x_grid});
                                [c_fbr,~] = hist3 ([state.mom{1}, state.mom{2}], 'Ctrs',{space.dof{1}.p_grid space.dof{2}.p_grid});
                            end
                            c_wig = 0; % dummy
                        otherwise
                            prt.error ('Binning of trajectories in more than 2 dimensions not yet implemented. Should try reduced_1d|2d instead')
                    end
                    
                    dvr_max = max ( dvr_max, max(c_dvr(:)) );
                    fbr_max = max ( fbr_max, max(c_fbr(:)) );
                    wig_max = max ( wig_max, max(c_wig(:)) );
                    
                end
            end
            
        else
            
            dvr_max = zeros (space.n_dim,1);
            fbr_max = zeros (space.n_dim,1);
            wig_max = zeros (space.n_dim,1);
            
            for m=1:hamilt.coupling.n_eqs
                if expect.pop.cha{m}(step)>expect.min_pop
                    
                    for k=1:space.n_dim
                        if strcmpi(info.system,'Matlab') % histcounts[2] not available in Octave
                            [c_dvr,~] = histcounts  (state.pos{1},               space.dof{1}.x_grid                     );
                            [c_fbr,~] = histcounts  (              state.mom{1},                      space.dof{1}.p_grid);
                            [c_wig,~] = histcounts2 (state.pos{1}, state.mom{1}, space.dof{1}.x_grid, space.dof{1}.p_grid);
                        elseif strcmpi(info.system,'Octave') % use of histc/hist3 is discouraged in Matlab
                            c_dvr = histc ( state.pos{1},                         space.dof{1}.x_grid                      );
                            c_fbr = histc (               state.mom{1},                                space.dof{1}.p_grid );
                            c_wig = hist3 ([state.pos{1}, state.mom{1}], 'Ctrs', {space.dof{1}.x_grid, space.dof{1}.p_grid});
                        end
                        dvr_max(k) = max ( dvr_max(k), max(c_dvr(:)) );
                        fbr_max(k) = max ( fbr_max(k), max(c_fbr(:)) );
                        wig_max(k) = max ( wig_max(k), max(c_wig(:)) );
                    end
                    
                end
            end
            
        end
        
        if isempty(plots.density.scale_dvr)
            plots.density.scale_dvr = dvr_max;
        end
        if isempty(plots.density.scale_fbr)
            plots.density.scale_fbr = fbr_max;
        end
        if isempty(plots.density.scale_wig)
            plots.density.scale_wig = wig_max;
        end
        
    end
    
    %% Q/M only: Determine ranges of kinetic/potential/total energy
    if isa(state,'wave') || ( isa(state,'traj') && space.n_dim < 3 )
        if isempty(plots.density.kin_min)
            plots.density.kin_min = 0;
        end
        if isempty(plots.density.kin_max)
            plots.density.kin_max = hamilt.kin_max;
        end
        plots.density.kin_delta = plots.density.kin_max - plots.density.kin_min;
        
        if isempty(plots.density.pot_min)
            plots.density.pot_min = hamilt.pot_min;
        end
        if isempty(plots.density.pot_max)
            plots.density.pot_max = hamilt.pot_max;
        end
        if plots.density.pot_min==plots.density.pot_max % Dirty trick for free particle
            plots.density.pot_max = plots.density.kin_max;
        end
        plots.density.pot_delta = plots.density.pot_max - plots.density.pot_min;
        
        plots.density.tef_min = min(plots.density.pot_min,0);
        plots.density.tef_max = plots.density.pot_max + plots.density.kin_max;
        plots.density.tef_delta = plots.density.tef_max - plots.density.tef_min;
    end
    
    % Determine maximal values of densities in pos/mom representation
    if isa (state,'wave') && space.n_dim <= 3
        
        dvr_max = 0;
        fbr_max = 0;
        
        for m=1:hamilt.coupling.n_eqs
            if expect.pop.cha{m}(step)>expect.min_pop
                dvr_max = max ( dvr_max, max(abs(state.dvr{m}(:)).^2) );
                
                fbr = state.dvr{m};
                for k = 1:space.n_dim
                    fbr = dvr2fbr(space.dof{k}, fbr);
                end
                fbr_max = max ( fbr_max, max(abs(fbr(:)).^2) );
                
            end
        end
        if isempty(plots.density.scale_dvr)
            plots.density.scale_dvr = dvr_max;
        end
        if isempty(plots.density.scale_fbr)
            plots.density.scale_fbr = fbr_max;
        end
        
    end
    
end

%% Animate densities in DVR/FBR/phase space

% Toggle plotting
if isfield (plots,'density')
    
    % First figure
    if step>1
        figure(1);
    
    % First call only
    else
        
        % Set figure size
        if isfield (plots,'expect')
            window_width = plots.density.w_width + plots.expect.w_width;
            plots.density.wide = true;
        else
            window_width = plots.density.w_width;
            plots.density.wide = false;
        end
        
        if strcmpi(info.system,'Matlab')
            h1 = figure(1);
            set(h1, ...
                'Name','Densities', ...
                'units','pixels', ...
                'position',[...
                plots.density.w_left ...
                plots.density.w_lower ...
                window_width ...
                plots.density.w_height] );
        elseif strcmpi(info.system,'Octave')
            figure(1,...
                'Name','Densities', ...
                'units','pixels', ...
                'position',[...
                plots.density.w_left ...
                plots.density.w_lower ...
                window_width ...
                plots.density.w_height] );
        end
        clf
        
        % Logos in all four corners of the plots
        if plots.density.logo
            show_logo (plots.density)
        end
        
        % Initialize movie export (first step only)
        % For the while being, this is not available in Octave
        if plots.density.export 
            
            if isempty(plots.density.file)
                state_type = state.string0;
                plot_type = class(plots.density);
                plot_type = plot_type (5:end);
                plots.density.file = strcat(state_type,'_',plot_type);
            end
            
            prt.disp ('***************************************************************')
            if ~plots.density.images && strcmpi(info.system,'Matlab')
                if ispc||ismac
                    extension = '.mp4';
                    myprofile = 'MPEG-4';
                elseif isunix
                    extension = '.avi';
                    myprofile = 'Motion JPEG AVI';
                end
                prt.disp ( ['Creating animated density plot : ' strcat(plots.density.file,extension)] )
                if exist (strcat(plots.density.file,extension), 'file')
                    delete( strcat(plots.density.file,extension) );
                end
                close(writerObj);
                writerObj = VideoWriter (plots.density.file, myprofile);
                open(writerObj);
            else
                prt.disp( 'Creating animated density plot as sequence of jpeg files')
            end
            prt.disp ('***************************************************************')
        end
    end
    
    % Various types of plots (with/without marginals)
    show_plot (plots.density, state, step)
    
    % Draw also expectation values (optionally)
    if isfield (plots,'expect')
        plots.expect.wide = true;
        show_expect (plots.expect,step);
        show_detail (plots.expect,state)
    end
    
    % Info about rendering, double buffering (not available in Octave)
    if step==1 && strcmpi(info.system,'Matlab')
        r = get (gcf, 'Renderer');
        d = get (gcf, 'DoubleBuffer');
        prt.disp (' ')
        prt.disp (['Type of density plot             : ' class(plots.density) ])
        prt.disp (['Rendering method                 : ' r])
        prt.disp (['Double buffering (painters only) : ' d])
        prt.disp (' ')
    end
    
    % Save last snapshot
    if step==time.steps.m_number
        if plots.density.export
            full_file_name1 = strcat(plots.density.file, '.jpg');
            full_file_name2 = strcat(plots.density.file, '.fig');
            prt.disp ( '***************************************************************' )
            prt.disp ( ['Saving last snapshot to file : ' full_file_name1] )
            prt.disp ( ['Saving last snapshot to file : ' full_file_name2] )
            prt.disp ( '***************************************************************' )
            prt.disp (' ')
            saveas(gcf,full_file_name1)
            saveas(gcf,full_file_name2)
        end
    end
    
    % Capture the current axes and add them as a movie frame. 
    % If we export single images, we output a new file instead
    if plots.density.export 
        if ~plots.density.images && strcmpi(info.system,'Matlab')
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        else
            % pad the filename with as many zeros as required
            padding = floor(log10(time.steps.m_number)) + 1;
            pattern = strcat('%s%0', num2str(padding), 'i.jpg');
            full_file_name = sprintf(pattern, plots.density.file, step);
            saveas(gcf, full_file_name);
        end
    end
    
    % Close/clear movie object (last step only)
    if plots.density.export && step==time.steps.m_number && strcmpi(info.system,'Matlab')
        if ~plots.density.images
            close(writerObj);
        end
    end
    
end

%% Animate expectation values (with uncertainties) if no densities are displayed

% Toggle plotting
if isfield(plots,'expect') && ~isfield(plots,'density')
    
    % Second figure
    if step>1
        figure(2);
    
    % First call only ...
    else
        
        % Clear figure and set size
        if strcmpi(info.system,'Matlab')
            h2 = figure(2);
            set(h2, ...
                'Name','Expectation values', ...
                'units','pixels', ...
                'position', [...
                plots.expect.w_left ...
                plots.expect.w_lower ...
                plots.expect.w_width ...
                plots.expect.w_height] );
        elseif strcmpi(info.system,'Octave')
            figure(2,...
                'Name','Expectation values', ...
                'units','pixels', ...
                'position', [...
                plots.expect.w_left ...
                plots.expect.w_lower ...
                plots.expect.w_width ...
                plots.expect.w_height] );
        end
        
        % Logos in the corners of the plots
        if plots.expect.logo
            show_logo (plots.expect)
        end
        
    end
    
    % Draw expectation values
    plots.expect.wide = false;
    show_expect (plots.expect);
    show_detail (plots.expect,state)
    
    % Last step only
    if step==time.steps.m_number
        
        % Export graphics to file (if desired)
        if plots.expect.export
            if isempty(plots.expect.file)
                state_type = state.string0;
                plot_type = class(plots.expect);
                plot_type = plot_type (5:end);
                plots.expect.file = strcat(state_type,'_',plot_type);
            end
            full_file_name1 = strcat(plots.expect.file, '.jpg');
            full_file_name2 = strcat(plots.expect.file, '.fig');
            prt.disp ( '***************************************************************' )
            prt.disp ( ['Saving plots of expectation values: ' full_file_name1] )
            prt.disp ( ['Saving plots of expectation values: ' full_file_name2] )
            prt.disp ( '***************************************************************' )
            prt.disp ( ' ' )
            saveas(gcf,full_file_name1)
            saveas(gcf,full_file_name2)
        end
        
    end
    
    drawnow;
end


%% Q/M only: Absorption spectrum 
if isa (state,'wave')
    
    % Toggle plotting (last step only, TDSE only)
    if isfield (plots,'spectrum') && step==time.steps.m_number && step>1
        
        % Fourier transform of autocorrelation
        spectrum (time.steps);
        
        % Third figure
        if strcmpi(info.system,'Matlab')
            h3 = figure(3);
            set(h3,...
                'Name','Spectrum', ...
            'units','pixels', ...
            'position', [...
                plots.spectrum.w_left ...
                plots.spectrum.w_lower ...
                plots.spectrum.w_width ...
                plots.spectrum.w_height] );
        elseif strcmpi(info.system,'Octave')
            figure(3,...
                'Name','Spectrum', ...
            'units','pixels', ...
            'position', [...
                plots.spectrum.w_left ...
                plots.spectrum.w_lower ...
                plots.spectrum.w_width ...
                plots.spectrum.w_height] );
        end
        
        % Logos in the corners of the plots
        if plots.spectrum.logo
            show_logo (plots.spectrum)
        end
        
        % Draw spectrum
        show_spec (plots.spectrum);
        
        % Export graphics to file (if desired)
        if plots.spectrum.export
            if isempty(plots.spectrum.file)
                state_type = state.string0;
                plot_type = class(plots.spectrum);
                plot_type = plot_type (5:end);
                plots.spectrum.file = strcat(state_type,'_',plot_type);
            end
            full_file_name1 = strcat(plots.spectrum.file, '.jpg');
            full_file_name2 = strcat(plots.spectrum.file, '.fig');
            prt.disp ( '***************************************************************' )
            prt.disp ( ['Saving plot of spectrum to file : ' full_file_name1] )
            prt.disp ( ['Saving plot of spectrum to file : ' full_file_name2] )
            prt.disp ( '***************************************************************' )
            prt.disp ( ' ' )
            saveas(gcf,full_file_name1)
            saveas(gcf,full_file_name2)
        end
        
        drawnow;
    end
end
