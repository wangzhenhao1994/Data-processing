
%--------------------------------------------------------------------------    
% Everything about expectation values and their uncertainties
%--------------------------------------------------------------------------
    
% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% see the README file for license details.
classdef generic < handle
    
    properties (Access=public)
        cha         % mean/average values for each channel
        unc         % uncertainties for each channel
        tot         % total, i.e. averaged over all channels
        qua         % which quantity (three letter code)
        ind         % indexing used for DVR, FBR or of AMOs
    end
    
    methods  (Access = public)
        
        % Constructor: Set default values
        function obj = generic (quantity)
            global hamilt time
            
            obj.cha = cell ( hamilt.coupling.n_eqs, 1);
            obj.unc = cell ( hamilt.coupling.n_eqs, 1);
            
            for m = 1:hamilt.coupling.n_eqs
                obj.cha{m} = zeros ( time.steps.m_number, 1);
                obj.unc{m} = zeros ( time.steps.m_number, 1);
            end
            
            obj.tot = 0;
            
            switch quantity
                case 'pop'
                case 'amo'
                case 'pos'
                case 'mom'
                case 'pot'
                case 'kin'
                case 'pur'
                otherwise
                    prt.error ('Invalid choice of quantity for mean value / uncertainty')
            end
            obj.qua = quantity;
            
        end
        
        traj_cha (obj, state, step)
        traj_tot (obj, state, step)
        
        wave_cha (obj, state, step)
        wave_fbr (obj, state, step)
        wave_kin (obj, state, step)
        wave_pop (obj, state, step)
        wave_tot (obj, state, step)
        
    end
    
end