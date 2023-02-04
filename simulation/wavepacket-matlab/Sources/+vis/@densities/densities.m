%------------------------------------------------------------------------------
%
% Settings for all density plots
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2012 Ulf Lorenz
%
% see the README file for license details.

classdef densities < vis.styles & handle
    
    properties (Access = public)
        
        w_left      % Left border of plot window
        w_lower     % Left border of plot window
        w_width     % Width of plot window
        w_height    % Height of plot window

        scale_dvr   % Scaling factor: position space 
        scale_fbr   % Scaling factor: momentum space 
        scale_wig   % Scaling factor: phase space
        
        hold        % Hold previous densities or overwrite them
        represent   % Position (dvr) or momentum (fbr) densities
        complex     % Choice of real|imag|abs2
        wide        % Toggle wide format: 16:9
        logo        % Toggle logos in all four corners
        energy      % Show also energy functions
        
        pot_min     % Minimum of potential energy
        pot_max     % Maximum of potential energy
        pot_delta   % Range of potential energy  
        
        kin_min     % Minimum of kinetic energy
        kin_max     % Maximum of kinetic energy
        kin_delta   % Range of kinetic energy
        
        tef_min     % Minimum of total energy 
        tef_max     % Maximum of total energy 
        tef_delta   % Range of total energy 
        
        range       % Manual setting of plot range
        x_min       % Minimum of first coordinate
        x_max       % Maximum of first coordinate
        y_min       % Minimum of second/momentum coordinate
        y_max       % Maximum of second/momentum coordinate
        z_min       % Minimum of third coordinate: height
        z_max       % Maximum of third coordinate: height
        
        export      % Toggle video export
        file        % Set the output to a custom filename (w/o suffix)
        images      % Export movie as a series of images
        
        mask        % Find steps where individual populations are non-negligible
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = densities
            
            obj = obj@vis.styles;
            
            % Position and size of plot window: Format 9:9
            obj.w_left   = round(  obj.s_height/16);   % Left border
            obj.w_lower  = round(  obj.s_height/16);   % Lower border
            obj.w_width  = round(9*obj.s_height/16);   % Width
            obj.w_height = round(9*obj.s_height/16);   % Height
            obj.w_width =  2*round(obj.w_width /2);    % Should be even
            obj.w_height = 2*round(obj.w_height/2);    % Should be even
  
            % Factor by which densities are divided 
            obj.scale_dvr = [];
            obj.scale_fbr = [];
            obj.scale_wig = [];
            
            % Various representation issues
            obj.hold = false;            % Hold previous densities or overwrite them
            obj.represent = 'dvr';       % Position (dvr) or momentum (fbr) densities
            obj.complex='abs2';          % Choice of real|imag|abs2
            obj.wide = [];               % Square or rectangular ("wide" 16:9) window
            obj.logo = true;             % Logos in all four corners
            obj.energy = true;           % Show also energy functions
            
            % Minimum and maximum ranges of kinetic/potential/total energy
            obj.kin_min = [];
            obj.kin_max = [];
            obj.kin_delta = [];
            
            % Minimum and maximum ranges of kinetic/potential energy
            obj.pot_min = [];
            obj.pot_max = [];
            obj.pot_delta = [];
            
            % Minimum and maximum ranges of total energy
            obj.tef_min = [];
            obj.tef_max = [];
            obj.tef_delta = [];
            
            % Manual setting of plot range
            obj.range = false;           % Toggle manual setting of plot range
            obj.x_min = [];              % Minimum of first coordinate
            obj.x_max = [];              % Maximum of first coordinate
            obj.y_min = [];              % Minimum of second/momentum coordinate
            obj.y_max = [];              % Maximum of second/momentum coordinate
            obj.z_min = [];              % Minimum of third coordinate: height
            obj.z_max = [];              % Maximum of third coordinate: height
            
            % Export to MP4 video file
            obj.export = true;           % Toggle video export
            obj.file   = [];             % Set the output to a custom filename (w/o suffix)
            obj.images = false;          % Export movie as a series of images
            
        end
        
        show_colors (obj) % see separate file
        
    end
    
end

