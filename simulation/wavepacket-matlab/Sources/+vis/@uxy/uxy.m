%--------------------------------------------------------------------------
%
% Plot evolution of input/state/output of bilinear control problem 
% in one figure with 2 or 3 subplots 
% (1) input=control field(s) u(t) - not drawn for the field-free case - 
% (2) state vector x(t)
% (3) output=observable(s) y(t) 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

classdef uxy < vis.styles & handle
    
    properties (Access = public)
        
        w_left      % Left border of plot window
        w_lower     % Left border of plot window
        w_width     % Width of plot window
        w_height    % Height of plot window
        
        legends     % Toggle curve legends
        
    end
    
    properties (Access = private)
        
        writerObj   % Movie file
        numPlots    % Two or three subplots
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = uxy
            
            obj = obj@vis.styles;
            
            % Position and size of plot window
            obj.w_left   = round(   obj.s_height/16);         % Left border
            obj.w_lower  = round(   obj.s_height/16);         % Lower border
            obj.w_width  = round( 7*obj.s_height/16);         % Width
            obj.w_height = round(12*obj.s_height/16);         % Height
            obj.w_width =  2*round(obj.w_width /2);           % Should be even
            obj.w_height = 2*round(obj.w_height/2);           % Should be even

            obj.legends = false;
        end
        
        initial  ( obj )  % open figure/subplots and prepare legends
        equilib  ( obj )  % draw equilibrium values (dashed lines)
        forward  ( obj,  first, last )  % draw forward propagation (solid lines)
        backward ( obj,  first, last )  % draw backward propagation (dotted lines)
        clearfig ( obj )  % clear figure/subplots
        closefig ( obj )  % close figure/subplots

    end
    
end
