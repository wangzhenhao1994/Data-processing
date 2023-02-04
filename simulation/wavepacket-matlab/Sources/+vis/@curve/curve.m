%--------------------------------------------------------------------------
%
% Visualize densities either in DVR ("coordinate" space) or FBR
% ("momentum" space). This plot accepts essentially all kinds of 
% grids, but note that it does not look good with all of them.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2009 Ulf Lorenz
%
% see the README file for license details.

classdef curve < vis.densities & handle
        
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = curve
            
            obj = obj@vis.densities;
            
        end
        
        %-------------------------------------------------
        % Show densities (from quantum wave functions)
        %-------------------------------------------------
        function show_plot (obj, state, step)
            
            global info space
            
            if ~isa(state,'wave') && ~isa(state,'traj')
                prt.error('Curve plots available only for wavefunctions or trajectory bundles')
            end
            
            if space.n_dim > 1
                prt.error('Cannot draw a curve plot for more than 1 dimension')
            end
            
            if obj.wide % Wide format: 16:9
                w=16; h=09;
            else % Square format: 9:9
                w=09; h=09;
            end
            
            if isa(state,'wave') && strcmpi(obj.complex,'abs2')
                if strcmpi (info.system,'Matlab')
                    subplot ( 'Position', [1/w 1/h 7/w 1/h] )
                    show_colors (obj)
                end
                subplot ( 'Position', [1/w 3/h 7/w 5/h] );
            else
                subplot ( 'Position', [1/w 1/h 7/w 7/h] );
            end
            
            if ~obj.hold
                hold off
            end
            plot( 1234567890, 1234567890 );
            hold on;
            
            switch lower(obj.represent)
                case 'dvr'
                    show_1d_dvr ( obj, state, step );
                case 'fbr'
                    show_1d_fbr ( obj, state, step );
                otherwise
                    prt.error ('Wrong choice of representation for curve plots')
            end
            
        end
        
    end
    
    methods (Access = private)
        
        show_1d_dvr ( obj, state, step )
        show_1d_fbr ( obj, state, step )

    end
    
end

