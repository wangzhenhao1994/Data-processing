%--------------------------------------------------------------------------
%
% Plot evolution of input/output of bilinear control problem 
% in one figure with 1 or 2 subplots 
% (1) input=control field(s) u(t) - not drawn for the field-free case - 
% (2) output=observable(s) y(t) 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-20.. Burkhard Schmidt
%
% see the README file for license details.

classdef in_out < vis.styles & handle
    
    properties (Access = public)
        
        w_left     % Left border of plot window
        w_lower    % Left border of plot window
        w_width    % Width of plot window
        w_height   % Height of plot window

        legends     % Toggle curve legends
        logo        % Toggle logos in all four corners
        wide        % Toggle wide format: 16:9
        
        u_min       % minimum value of control field(s) u(t)
        u_max       % maximum value of control field(s) u(t)
        
        y_min       % minimum value of observable(s) u(t)
        y_max       % maximum value of observable(s) u(t)
        
        export      % Toggle graphics export
        file        % Set the output to a custom filename (w/o suffix)
        
    end
    
    properties (Access = private)
        
        
    end

    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = in_out
            
            % Inheriting from superclass
            obj = obj@vis.styles;
            
            % Position and size of plot window: Format 7:9
            obj.w_left   = round(11*obj.s_height/16);  % Left border
            obj.w_lower  = round(01*obj.s_height/16);  % Lower border
            obj.w_width  = round(07*obj.s_height/16);  % Width
            obj.w_height = round(09*obj.s_height/16);  % Height
            obj.w_width =  2*round(obj.w_width /2);    % Should be even
            obj.w_height = 2*round(obj.w_height/2);    % Should be even

            % Appearence of plot
            obj.legends  = true;         % Toggle legends
            obj.logo     = true;         % Toggle logos in all four corners
            obj.wide = [];               % Toggle wide format: 16:9

            % Range settings for input fields
            obj.u_min = -1;
            obj.u_max = +1;
            
            % Range settings for output observables
            obj.y_min = -1;
            obj.y_max = +1;
            
            
        end
        
        %-------------------------------------------------
        % Show densities (from quantum wave functions)
        %-------------------------------------------------
        function show_expect (obj,step)
            
            global info state time
            
            switch lower(class(state))
                case {'ket','rho'}
                otherwise
                    prt.error('"in_out" plots available only for ket vectors or rho matrices')
            end            
            
            if obj.wide % Wide format: 16:9
                w=16; h=09; o=9/16;
            else % Narrow format: 7:9
                w=07; h=09; o=0;
            end
            
            % Two subplots
            if strcmpi (info.program,'qm_propa') && isfield(time,'pulse')
                
                % Input (control field(s) versus time
                subplot ( 'Position', [o+1/w 5/h 5/w 3/h] )
                cla;
                input (obj,step)
                if obj.wide
                    set ( gca, 'YAxisLocation', 'right');
                end
                
                % Output: observables versus time
                subplot ( 'Position', [o+1/w 1/h 5/w 3/h] )
                cla;
                output (obj,step)
                if obj.wide
                    set ( gca, 'YAxisLocation', 'right');
                end
                
                % One subplot
            else
                
                % Output: observables versus time
                subplot ( 'Position', [o+1/w 1/h 5/w 7/h] )
                output (obj,step)
                if obj.wide
                    set ( gca, 'YAxisLocation', 'right');
                end
                
            end
        end
        
        % Currently, this is just a dummy
        function show_detail (~,~)
        end
        
    end
    
    methods (Access = private)
        input(obj,step)
        output(obj,step)
    end
    
end

