%--------------------------------------------------------------------------
% Dummy class for empty system bath coupling
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2020-.... Burkhard Schmidt's group
%
% see the README file for license details.

% The only reason why we need this class is because we want to 
% be able to identify vanishing system bath coupling
% unequivocally by isa(hamilt.sbc{m,n},'sbc.empty')
classdef empty < sbc.generic
            
    methods (Access = public)
        
        function disp (obj)
            disp@sbc.generic(obj)
            prt.disp ('Not available')
            prt.disp ('***************************************************************')
            prt.disp (' ')
        end
        
    end
    
end

