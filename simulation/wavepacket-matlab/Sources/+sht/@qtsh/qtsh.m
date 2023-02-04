%--------------------------------------------------------------------------
%
% Quantum trajectory surface hopping
% by Craig C. Martens
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2020-.... Leonardo Cancissu Araujo
%
% see the README file for license details.

classdef qtsh < sht.fssh & handle
    
    methods (Access = public)

        % Constructor: Setting defaults and a text string
        function obj = qtsh (n,seed)
            obj = obj@sht.fssh(n,seed);     % Inherit from superclass
            obj.string0 = 'qtsh';
            obj.string4 = 'Quantum trajectory surface hopping';
            
            obj.choice_eq_motion = 'qt_adi';    % quantum trajectories: adiabatic
            obj.rescale  = false;               % do not rescale momentum after a hop
        end   
        
    end
end

