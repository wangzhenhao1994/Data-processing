%--------------------------------------------------------------------------
% Dummy class for empty polarizability
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2020-.... Burkhard Schmidt's group
%
% see the README file for license details.

% The only reason why we need this class is because we want to 
% be able to identify vanishing polarizability
% unequivocally by isa(hamilt.pol{m,n},'pol.empty')
classdef empty < pol.generic
            
    methods (Access = public)
        
        function disp (obj)
            disp@pol.generic(obj)
            prt.disp ('Not available')
            prt.disp ('***************************************************************')
            prt.disp (' ')
        end
        
    end
    
end

