%--------------------------------------------------------------------------
%
% MSSH = Multiple switches surface hopping (historic original)
%
% Determine transition probability from  |c(t)|^2
% from TDSEs attached to each of the trajectories.
% While in principle correct, this SH variant is 
% known to switch far too often, even when the
% trajectories are outside the regions of strong 
% non-adiabatic coupling
%
% see: J. C. Tully, R. K. Preston
%      J. Chem. Phys. 55(2), 562-572 (1971)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef mssh < sht.generic & handle
    
    methods (Access = public)
        
        % Constructor: Setting defaults and a text string
        function obj = mssh (n,seed)
            obj = obj@sht.generic(n,seed);     % Inherit from superclass
            obj.string0 = 'mssh';
            obj.string4 = 'Multiple switches surface hopping';
            
            obj.extra_tdse_solving = 1; % Solve the TDSE if it is not already done by the choice of propagator
        end   
        
            
        % Get probabilities of hopping from state "m" to state "n"
        function probable = prob_hop (obj,m,n,ind_m)
            
            % Multiple switches surface hopping (historic original)
            probable = abs(obj.psi{n}(ind_m)).^2;
            
        end

        
    end
    
end

