% Kinetic energy in Jacobi coordinates for the angle.
% See for example J. Chem. Phys 116:4403
 
% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%               2017 Burkhard Schmidt
%
% see the README file for license details.

classdef jacobi < handle
       
    properties (Access = public)
        
        r_0         % Constant distance ...
        dof_R       % index of dof: ... distance ...
        dof_r       % index of dof: ... distance ...
        dof_c       % index of dof: ... angle ...
        mass_R      % Reduced mass of ...
        mass_r      % Reduced mass of ...
        
    end
    
    properties (Access = private)
        
        grid        % DVR grid representation of this kinetic energy operator
        grid_exp    % same as grid, but exponentiated, for split operator
        
    end
    
    methods (Access = public)
        
        kinetic     (obj, psi, new)
        kinetic_exp (obj, psi)
        init_kin    (obj, fraction, output)
        dvrkin = kinetic2dvr(obj)
        obj = subsasgn(obj, index, val)
        retval = subsref(obj, s)
        
    end
    
end
