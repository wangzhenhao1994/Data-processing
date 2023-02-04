%--------------------------------------------------------------------------
% Dummy class for empty potential energy function: Free particle(s)
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2020-.... Burkhard Schmidt's group
%
% see the README file for license details.

% The only reason why we need this class is because we want to 
% be able to identify vanishing potential energy functions
% unequivocally by isa(hamilt.pot{m,n},'pot.empty')
classdef empty < pot.generic
    
    methods (Access = public)
        
        % Display potential, overloading default disp method
        function disp (obj)
            disp@pot.generic (obj)
            prt.disp ('Not available')
            prt.disp ('***************************************************************')
            prt.disp (' ')
        end
        
    end
    
end

