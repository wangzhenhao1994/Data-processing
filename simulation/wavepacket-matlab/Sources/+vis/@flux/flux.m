%--------------------------------------------------------------------------
%
% Calculate and visualize flux density (1D)
%
% J = -i*hbar/2 * [       psi^* * grad (psi) - psi * grad (psi^*) ]
% 
%   =    hbar * Im [      psi^* * grad (psi) ]
%
%   =    hbar * Re [ -i * psi^* * grad (psi) ]
%
%   =           Re [      psi^* *   P  (psi) ]
%
% see, e.g., Cohen-Tannoudji, Chapter III, Complement D, Eq. D-17 (p.239)
%
% normalization with 1/mass ?!?
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2013-2017 Burkhard Schmidt
%
% see the README file for license details.

classdef flux < vis.densities & handle
    
    properties (Access = public)
        
        cnt_nlev    % Number of contours: density/energy
        cnt_levels  % Toggle manual settings of color ranges
        cnt_min     % Lower contour line for density plots
        cnt_max     % Upper contour line for density plots
     
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = flux
            
            obj = obj@vis.densities;
            
            obj.cnt_nlev  = [30 15];     % Number of contours: density/energy
            obj.cnt_levels = false;      % Toggle manual settings of color ranges
            obj.cnt_min = [];            % Lower contour line for density plots
            obj.cnt_max = [];            % Upper contour line for density plots

        end
 


        %---------------------------------------------------------
        % Show flux of densities 
        %---------------------------------------------------------
        function show_plot (obj, state, step)
            
            global space
            
            if ~isa(state,'wave') && ~isa(state,'traj')
                prt.error('Flux plots available only for wavefunctions or trajectory bundles')
            end

            if space.n_dim > 2
                prt.error('Cannot draw a flux plot for more than 2 dimensions')
            end
            
            if obj.wide % Wide format: 16:9
                w=16; h=09;
            else % Square format: 9:9
                w=09; h=09;
            end
            
            switch space.n_dim
                case  1
                    
                    % Lower part: flux density
                    subplot ( 'Position', [1/w 1/h 7/w 3/h] );
                    hold off; plot( 1234567890, 1234567890 ); hold on;
                    show_1d_lower ( obj, state, step )
                    
                    % Upper part: probability density
                    subplot ( 'Position', [1/w 5/h 7/w 3/h] );
                    hold off; plot( 1234567890, 1234567890 ); hold on;
                    show_1d_upper ( obj, state, step )
                    
                case  2
                    subplot ( 'Position', [1/w 1/h 7/w 7/h] );
                    show_2d ( obj, state, step )
                otherwise
                    prt.error('Cannot draw a flux plot for more than 2 dimensions')
            end
            
        end
        
    end
    
    methods (Access = private)
       
        show_1d_lower ( obj, state, step )
        show_1d_upper ( obj, state, step )
        show_2d       ( obj, state, step )
        
    end
    
end



