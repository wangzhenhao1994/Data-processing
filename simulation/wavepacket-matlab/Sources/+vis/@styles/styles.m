%------------------------------------------------------------------------------
%
% Style settings for all WavePacket plots
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2012 Ulf Lorenz
%
% see the README file for license details.

classdef styles < handle
    
    properties (Access = public)
        
        colors      % RGB colors for curves/lines
        patterns    % Patterns for curves/lines
        
        l_thin      % Widths for thin curves/lines
        l_thick     % Widths for thick curves/lines
        l_extra     % Widths for extrathick curves/lines
        
        f_small     % Size of small fonts
        f_large     % Size of large fonts
        f_heavy     % Weight of heavy fonts
        f_light     % Weight of light fonts
        f_name      % Font name
        
        m_small     % Size of small markers
        m_large     % Size of large markers
        m_extra     % Size of extra large markers
        m_huger     % Size of even huger markers
        
        s_height    % Screen height
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = styles
            
            screen = get (0, 'ScreenSize');  % Screen size
            obj.s_height = screen(4);    % Don't use width for double screens
            
            % If matlab is started from the command line (matlab -nodisplay), the
            % above code screws up. In this case, supply a standard value for the
            % dimensions.
            if obj.s_height < 100
                obj.s_height = 768;
            end
            
            % RGB colors for individual channels: Matlab's color scheme introduced in R2014a
            obj.colors = [...
                0    0.4470    0.7410; ...
                0.8500    0.3250    0.0980; ...
                0.9290    0.6940    0.1250; ...
                0.4940    0.1840    0.5560; ...
                0.4660    0.6740    0.1880; ...
                0.3010    0.7450    0.9330; ...
                0.6350    0.0780    0.1840];
            obj.colors = repmat(obj.colors,10,1); % repeat 10 times
            
            % Patterns for curves/lines
            obj.patterns = {...
                '-'; ...                 % 1: solid (default)
                '--';  ...               % 2: dashed
                ':'; ...                 % 3: dotted
                '-.'; ...                % 4: dash-dotted
                '-'; ...                 % 5: solid (default)
                '--';  ...               % 6: dashed
                ':'; ...                 % 7: dotted
                '-.'; ...                % 8: dash-dotted
                '-'; ...                 % 9: solid (default)
                '--';  ...               % 10: dashed
                ':'; ...                 % 11: dotted
                '-.'};                   % 12: dash-dotted
            
            % Line widths
            obj.l_thin  = 1;          % energy contours/curves
            obj.l_thick = 2;          % axes, density contours/curves
            obj.l_extra = 4;          % extrathick, for polar plots
            
            % Font sizes and weights
            obj.f_large = round(obj.s_height/125);
            obj.f_small = round(obj.s_height/125)-2;
            obj.f_heavy = 'bold';
            obj.f_light = 'light';
            obj.f_name  = 'Helvetica';
            
            % Marker sizes
            obj.m_small = round(obj.s_height/400)-2;
            obj.m_large = round(obj.s_height/400);
            obj.m_extra = round(obj.s_height/050);
            obj.m_huger = round(obj.s_height/025);
            
        end
        
        show_logo (obj) % see separate file
        
    end
    
    methods (Static)
        
        plot_color(x,y,c,c0,w,offset,swap) % see separate file
        
    end
    
end