%--------------------------------------------------------------------------
%
% Generic properties of all class definitions for initial wavefunctions
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.


classdef generic < handle
    
    properties (Access = public) 
        
        dof         % Index of degree of freedom
        dvr         % grid representation of this wavefunction

    end
    
end

