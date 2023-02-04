%--------------------------------------------------------------------------
%
% Visualize multidimensional wavepacket or trajectory propagations
% by means of one-dimensional reduced densities which are plotted 
% either in DVR or FBR or Wigner or R-R' representation
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-20xy Burkhard Schmidt's group
%               2007-2010 Ulf Lorenz
%
% see the README file for license details.

classdef reduced_1d < vis.densities & handle
    
    properties (Access = public)

        cnt_nlev    % Number of contours: density
        rho_max;    % Maximal values of (coupled) densities 
        
        expect      % Plot "trajectory" of mean values
        
    end
                
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = reduced_1d

            obj = obj@vis.densities;
            
            obj.cnt_nlev  = 20;          % Number of contours for reduced densities
            obj.rho_max;                 % Maximal values of (coupled) densities 
            obj.expect = false;          % Plot "trajectory" of mean values
            
        end
        
        %-------------------------------------------------
        % Show reduced densities
        %-------------------------------------------------
        function show_plot (obj, state, step)
            
            global expect space hamilt
            
            if ~isa(state,'wave') && ~isa(state,'traj')
                prt.error('Reduced density (1D) plots available only for wavefunctions or trajectory bundles')
            end
            
            %% Full Q/M only
            if isa(state,'wave')
                
                % Initialize calculations of purity
                for k = 1:space.n_dim
                    if step==1
                        expect.pur{k} = obs.generic ('pur');
                    end
                    expect.pur{k}.tot(step) = 0;
                end
                
                % Calculate reduced densities
                psi.redu = cell (hamilt.coupling.n_eqs, space.n_dim);
                
                % Loop over (coupled) wavefunctions
                for m = 1:hamilt.coupling.n_eqs
                    
                    % Only if population exceeds a certain threshold
                    if expect.pop.cha{m}(step)>expect.min_pop
                        
                        for k = 1:space.n_dim
                            psi.redu{m,k} = zeros(space.dof{k}.n_pts);
                        end
                        
                        switch space.n_dim
                            case 1
                                wave_redu_1d(obj,state,m);
                            case 2
                                wave_redu_2d(obj,state,m);
                            case 3
                                wave_redu_3d(obj,state,m);
                            case 4
                                wave_redu_4d(obj,state,m);
                            case 5
                                wave_redu_5d(obj,state,m);
                            case 6
                                wave_redu_6d(obj,state,m);
                            otherwise
                                prt.error ('Reduced density matrices only up to 6 dimensions')
                        end
                        
                        % Formally not really correct: As Max B. keeps telling us,
                        % the reduced densities should be normalized with the
                        % weights in N-1 dimensions. However, this error is comp-
                        % ensated in the calculation of the purities!
                        for k = 1:space.n_dim
                            psi.redu{m,k} = psi.redu{m,k} * space.weight(1);
                        end
                        
                        % Calculate purity as tr(rho^2)
                        expect.pur{k}.cha{m}(step) = abs(trace(psi.redu{m,k}^2));
                        expect.pur{k}.tot   (step) = expect.pur{k}.tot(step) ...
                            + expect.pur{k}.cha{m}(step) / expect.pop.cha{m}(step)^2;

                    end
                    
                end
                
            end
            
            %% Visualize reduced densities
            p = 1;
            m = ceil(sqrt(space.n_dim));
            for k = 1:space.n_dim
                
                % Arranging the subplots
                if obj.wide % Wide format: 16:9
                    subplot (m,2*m,p)
                    if mod(p,m)==0
                        p = p+1+m;
                    else
                        p = p+1;
                    end
                else % Square format: 9:9
                    subplot (m,m,p)
                    p = p+1;
                end
                if ~obj.hold; hold off; end
                plot( 1234567890, 1234567890 );
                hold on;
                
                % Making the individual plots
                switch(lower(obj.represent))
                    case 'dvr'
                        show_1d_dvr ( obj, state, step, k );
                    case 'fbr' % Q/C only
                        show_1d_fbr ( obj, state, step, k );
                    case 'wig'
                        show_1d_wig ( obj, state, step, k );
                    case 'rrp' % Q/M only
                        show_1d_rrp ( obj, state, step, k );
                    otherwise
                        prt.error ('Choose either one of DVR or FBR or WIG or RRP representation')
                end
                
            end
            
        end
    end
    
    methods ( Access = private )
        
        show_1d_dvr ( obj, state, step, k );
        show_1d_fbr ( obj, state, step, k );
        show_1d_wig ( obj, state, step, k );
        show_1d_rrp ( obj, state, step, k );
        
        wave_redu_1d ( obj, state, m)
        wave_redu_2d ( obj, state, m)
        wave_redu_3d ( obj, state, m)
        wave_redu_4d ( obj, state, m)
        wave_redu_5d ( obj, state, m)
        wave_redu_6d ( obj, state, m)
        
    end
end

