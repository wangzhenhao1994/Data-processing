% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2017 Burkhard Schmidt
%
% see the README file for license details.

classdef triatomic  < handle
    
    properties (Access = public)
        
        dof         % Indices of dof's for distances AB, BC
        mass        % Masses of atoms A, B, C
        theta       % Value of ABC angle 
        
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
