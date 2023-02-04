%--------------------------------------------------------------------------
%
% Visualize "ket" vectors as 2-D bar graph
% Visualize "rho" matrices as 3-D bar graph
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-.... Burkhard Schmidt
%
% see the README file for license details.

classdef bar < vis.densities & handle
        
    properties (Access = public)

        bar_view    % View point for surface plot: [az el]
        
    end

    
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = bar
            
            obj = obj@vis.densities;
            obj.bar_view  = [60 40];     % View point for surface plot: [az el]
            
        end
        
        %-------------------------------------------------
        % Show densities (from quantum wave functions)
        %-------------------------------------------------
        function show_plot (obj, state, step)
            
            switch lower(class(state))
                case {'ket','rho'}
                otherwise
                    prt.error('Bar plots available only for "ket" vectors or "rho" matrices')
            end
                       
            if obj.wide % Wide format: 16:9
                w=16; h=09;
            else % Square format: 9:9
                w=09; h=09;
            end
            
            subplot ( 'Position', [1/w 1/h 7/w 7/h] );
            
            if ~obj.hold
                hold off
            end
            plot( 1234567890, 1234567890 );
            hold on;
            
            switch lower(class(state))
                case 'ket'
                    show_ket ( obj, state, step );
                case 'rho'
                    if isempty (state.save_suffix)  % no balancing, truncation, etc
                        show_rho ( obj, state, step );
                    else
                        show_ket ( obj, state, step ); % with balancing, truncation, etc
                    end
            end
            
        end
        
    end
    
    methods (Access = private)
        
        show_ket ( obj, state, step )
        show_rho ( obj, state, step )

    end

    
end

