%--------------------------------------------------------------------------
%
% Polar plots of one-dimensional wavefunctions/densities on a periodic
% domain. The coordinate is the angle or its cosine, depending on the grid),
% the distance from the origin is the squared norm, the color is the phase.
% Since the plotting is a bit difficult to interpret, no potentials are
% ever plotted.
%
% Matlab offers a polar plot which we do not use for a couple of reasons.
% First, it does not really set up a polar plot, but just draws an "image"
% of a polar plot in an underlying x/y plot. You cannot really resize it etc.
% Second, it cannot be configured; zero angle is in the horizontal, which is
% very unintuitive at least for me. So it boils down to drawing a couple of
% lines that we can also do by hand.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2009 Ulf Lorenz
%               2009 Burkhard Schmidt
%
% see the README file for license details.

classdef polar < vis.densities & handle
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = polar
            
            obj = obj@vis.densities;
            
        end
        
        %-------------------------------------------------
        % Show densities (from quantum wave functions)
        %-------------------------------------------------
        function show_plot (obj, state, step)
            
            global space
            
            if ~isa(state,'wave')
                prt.error('Polar plots available only for wavefunctions')
            end

            if space.n_dim ~= 1
                prt.error('Cannot draw a polar plot for >1 dimensions.')
            end
            
            if obj.wide % Wide format: 16:9
                w=16; h=09;
            else % Square format: 9:9
                w=09; h=09;
            end
            
            subplot ( 'Position', [1/w 1/h 7/w 7/h] );
            hold off; plot( 1234567890, 123456789); hold on;
            
            switch lower(obj.represent)
                case 'dvr'
                    show_1d_dvr ( obj, state, step );
                otherwise
                    prt.error ('Wrong choice of representation for polar plots')
            end
        end
    end
    
    methods (Access = private)
           
        show_1d_dvr ( obj, state, step )
        
    end
    
end