%--------------------------------------------------------------
%
% Visualize 1-dim wavepacket (Wigner transform in phase space)
% with marginal distributions (i.e. position/momentum densities)
% or 2-dim wavepacket (either in DVR or in FBR representation)
% or 3-dim wavepacket (either in DVR or in FBR representation)
% using surface plots, with or without showing energy surfaces
%
%--------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

classdef surface < vis.densities & handle
    
    properties (Access = public)
        
        marginals   % Toggle marginal functions
        
        srf_view    % View point for surface plot: [az el]
        srf_look    % Look of surface plot [shading lighting]
        srf_light   % Angle of light for surface plot: [az el]
        srf_color   % Toggle manual settings of color ranges
        
        col_min     % Minimum energy value passed to caxis()
        col_max     % Maximum energy value passed to caxis()
        col_map     % Choice of Matlab's color maps
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = surface
            
            obj = obj@vis.densities;
            
            obj.marginals = true;        % Toggle marginal functions
            
            obj.srf_view  = [60 75];     % View point for surface plot: [az el]
            obj.srf_look  = [true true]; % Look of surface plot [shading lighting]
            obj.srf_light = [45 45];     % Angle of light for surface plot: [az el]
            obj.srf_color = false;       % Toggle manual settings of color ranges
            obj.col_min = [];            % Minimum energy value passed to caxis()
            obj.col_max = [];            % Maximum energy value passed to caxis()
            obj.col_map = 'default';     % Choice of one of Matlab's color maps
            
        end
        
        %-------------------------------------------------
        % Show densities 
        %-------------------------------------------------
        function show_plot (obj, state, step)
            
            global space
            
            if ~isa(state,'wave') && ~isa(state,'traj')
                prt.error('Surface plots available only for wavefunctions or trajectory bundles')
            end

            if space.n_dim > 3
                prt.error('Cannot draw a surface plot for more than 3 dimensions')
            end

            
            if obj.wide % Wide format: 16:9
                w=16; h=09;
            else % Square format: 9:9
                w=09; h=09;
            end
            subplot ( 'Position', [1/w 1/h 7/w 7/h] );
            
            switch space.n_dim
                case 1
                    vis.wigner  (      state, step )
                    show_1d_wig ( obj, state, step )
                case 2
                    if strcmpi(obj.represent,'dvr')
                        show_2d_dvr ( obj, state, step )
                    elseif strcmpi(obj.represent,'fbr')
                        show_2d_fbr ( obj, state, step )
                    end
                case 3
                    if isa (state,'traj')
                        prt.error('Surface plots for trajectory data in 3-dim not yet implemented')
                    end
                    if strcmpi(obj.represent,'dvr')
                        show_3d_dvr ( obj, state, step )
                    elseif strcmpi(obj.represent,'fbr')
                        show_3d_fbr ( obj, state, step )
                    end
                otherwise
                    prt.error('Cannot draw a surface plot for more than 3 dimensions')
            end
            
        end
        
    end
    
    methods (Access = private)
        
        show_1d_wig ( obj, state, step )
        show_2d_dvr ( obj, state, step )
        show_2d_fbr ( obj, state, step )
        show_3d_dvr ( obj, state, step )
        show_3d_fbr ( obj, state, step )
    end
    
end
            
            