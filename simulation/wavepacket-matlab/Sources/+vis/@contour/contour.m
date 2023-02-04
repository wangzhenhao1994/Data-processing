%--------------------------------------------------------------------------
%
% Visualize Wigner transform of 1-dim wavepacket in phase space 
% or classical phase space density from trajectory bundles
% using contour plots with 
% Optionally with marginals representing position and momentum densities. 
%
% Compose animated figure from three subplots (see included subfunctions)
% containing densities (position, momentum, phase space) and respective
% energy functions (potential, kinetic, total)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

classdef contour < vis.densities & handle
    
    properties (Access = public)
        
        marginals   % Toggle marginal function
        
        cnt_nlev    % Number of contours: density/energy
        cnt_levels  % Toggle manual settings of contour levels
        cnt_min     % Lower contour line for density plots
        cnt_max     % Upper contour line for density plots
        
        expect      % Plot "trajectory" of mean values
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = contour
            
            obj = obj@vis.densities;
            
            obj.marginals = true;        % Toggle marginal functions
            
            obj.cnt_nlev  = [30 15];     % Number of contours: density|energy
            obj.cnt_levels = false;      % Toggle manual settings of contour levels
            obj.cnt_min = [];            % Lower contour line for density plots
            obj.cnt_max = [];            % Upper contour line for density plots
            
            obj.expect = false;          % Plot "trajectory" of mean values
            
        end
        
        %-----------------------------------------------------
        % Show densities 
        %-----------------------------------------------------
        function show_plot (obj, state, step)
            
            global info space
            
            if ~isa(state,'wave') && ~isa(state,'traj')
                prt.error('Contour plots available only for wavefunctions or trajectory bundles')
            end

            if space.n_dim > 2
                prt.error('Cannot draw a contour plot for more than 2 dimensions')
            end
            
            if obj.wide % Wide format: 16:9
                w=16; h=09;
            else % Square format: 9:9
                w=09; h=09;
            end
            
            % Create Wigner distributions, optionally total energy function
            if space.n_dim == 1
                vis.wigner (state, step)
            end
            
            % Without marginals
            if ~obj.marginals
                
                % Wigner quasi-density, "mean trajectory" and total energy
                subplot ( 'Position', [1/w 1/h 7/w 7/h] )
                if ~obj.hold; hold off; end
                plot( 1234567890, 1234567890 );
                hold on;
                switch space.n_dim
                    case 1
                        show_1d_wig ( obj, state, step )
                    case 2
                        if strcmpi(obj.represent,'dvr')
                            show_2d_dvr ( obj, state, step )
                        elseif strcmpi(obj.represent,'fbr')
                            show_2d_fbr ( obj, state, step )
                        end
                end
                
                % With marginals
            else
                % Upper left: Wigner quasi-density, "mean trajectory" and total energy
                subplot ( 'Position', [1/w 3/h 5/w 5/h] )
                if ~obj.hold; hold off; end
                plot( 1234567890, 1234567890 );
                hold on;
                switch space.n_dim
                    case 1
                        show_1d_wig ( obj, state, step )
                    case 2
                        if strcmpi(obj.represent, 'dvr')
                            show_2d_dvr ( obj, state, step )
                        elseif strcmpi(obj.represent, 'fbr')
                            show_2d_fbr ( obj, state, step )
                        end
                end
                
                % Upper right: Momentum density and kinetic energy curve
                subplot ( 'Position', [6/w 3/h 2/w 5/h] )
                if ~obj.hold; hold off; end
                plot( 1234567890, 1234567890 );
                hold on;
                switch space.n_dim
                    case 1
                        show_1d_fbr  ( obj, state, step )
                    case 2
                        if strcmpi(obj.represent, 'dvr')
                            show_2d_dvr2  ( obj, state, step )
                        elseif strcmpi(obj.represent, 'fbr')
                            show_2d_fbr2  ( obj, state, step )
                        end
                end
                
                % Lower left: Position density and potential energy curve
                subplot ( 'Position', [1/w 1/h 5/w 2/h] )
                if ~obj.hold; hold off; end
                plot( 1234567890, 1234567890 );
                hold on;
                switch space.n_dim
                    case 1
                        show_1d_dvr ( obj, state, step )
                    case 2
                        if strcmpi(obj.represent, 'dvr')
                            show_2d_dvr1  ( obj, state, step )
                        elseif strcmpi(obj.represent, 'fbr')
                            show_2d_fbr1  ( obj, state, step )
                        end
                        
                end
                        
                % Lower right: Color mapping
                if isa (state,'wave') && space.n_dim == 1 && strcmpi(info.system,'Matlab')
                    subplot ( 'Position', [6/w 1/h 2/w 2/h] )
                    if ~obj.hold; hold off; end
                    plot( 1234567890, 1234567890 );
                    hold on;
                    if strcmpi (info.system,'Matlab')
                        show_colors ( obj )
                    end
                end
                
                
            end
            if strcmpi (info.system,'Octave')
                drawnow
            end
        end
        
    end
    
    methods (Access = private)

        show_1d_wig  ( obj, state, step )
        show_1d_dvr  ( obj, state, step )
        show_1d_fbr  ( obj, state, step )
        show_2d_dvr  ( obj, state, step )
        show_2d_dvr1 ( obj, state, step )
        show_2d_dvr2 ( obj, state, step )
        show_2d_fbr  ( obj, state, step )
        show_2d_fbr1 ( obj, state, step )
        show_2d_fbr2 ( obj, state, step )

    end
    
end

