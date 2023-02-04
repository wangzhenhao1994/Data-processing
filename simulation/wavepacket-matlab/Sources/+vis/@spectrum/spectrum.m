%----------------------------------------------------------
%
% Plot power spectrum (Fourier transform of autocorrelation)
%
%----------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%
% see the README file for license details.

classdef spectrum < vis.styles & handle
    
    properties (Access = public)
        
        w_left      % Left border of plot window
        w_lower     % Lower border of plot window
        w_width     % Width of plot window
        w_height    % Height of plot window
        
        logo        % Toggle logos in all four corners
              
        export      % Toggle graphics export
        file        % Set the output to a custom filename (w/o suffix)
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = spectrum
            
            obj = obj@vis.styles;
            
            % Position and size of plot window
            obj.w_left   = round(01  *obj.s_height/16);     % Left border
            obj.w_lower  = round(11.5*obj.s_height/16);     % Lower border
            obj.w_width  = round(16  *obj.s_height/16);     % Width
            obj.w_height = round(03  *obj.s_height/16);     % Height
            
            % Appearence of plot
            obj.logo     = true;         % Toggle logos in all four corners
            
            % Export to graphics file: jpg/tiff/epsc/pdf, etc ...
            obj.export   = false;        % Toggle graphics export
            obj.file = [];               % Set custom filename (suffix determines image type)
            
        end
        
        
        function show_spec (obj)
            
            global hamilt info time
            
            % Plot total population (black curve)
            if strcmp (info.program,'qm_bound')
                x = hamilt.eigen.eig_vals;
                y = ones(size(x));
                h = stem ( x,y );
                set(h, 'Marker', 'none')
            elseif strcmp (info.program,'qm_propa')
                h = plot ( time.steps.freq, abs(time.steps.spec));
            end
            set(h, ...
                'LineStyle', obj.patterns{1}, ...
                'Color',     'black', ...
                'LineWidth', obj.l_thick)
            
            % Axes range
            axis ( [ hamilt.pot_min hamilt.pot_max+hamilt.kin_max -0.1 1.1 ] )
            
            % Labels
            set(gca, ...
                'LineWidth',  obj.l_thick, ...
                'FontName',   obj.f_name, ...
                'FontSize',   obj.f_large, ...
                'FontWeight', obj.f_heavy)
            
            xlabel ( 'energy' )
            if info.program=='qm_propa'
                ylabel ( 'intensity' )
            end
            
        end
        
    end
    
end

