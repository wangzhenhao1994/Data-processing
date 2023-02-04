%****************************************************************
% Model function for transition dipole moment of H2O
% See J. Chem. Phys. 97:8285                        
%****************************************************************

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2009 Ulf Lorenz
%
% see the README file for license details.

classdef H2O < dip.generic & handle
    
    properties (Access = public)
        
        b           % range parameter
        r0          % equilibrium distance H-O
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = H2O
            obj.b  = 1;     % range parameter
            obj.r0 = 1.81;  % equilibrium distance H-O
        end
        
        % Initialize dipole moment: Set/check parameters
        function init (~)
            global hamilt space
            
            if space.n_dim ~= 2
                prt.error ('This dipole function is only for 2 dimensions')
            end
            
            if hamilt.coupling.n_eqs ~= 2
                prt.error ('This dipole function is only for 2 states')
            end
            
        end
        
        % Display dipole moment, overloading default disp method
        function disp(obj)
            disp@dip.generic(obj)
            prt.disp ('Model function for transition dipole moment of H2O')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp (['Range parameter          : ' num2str(obj.b)])
            prt.disp (['Equilibrium distance H-O : ' num2str(obj.r0)])
            prt.disp (' ')
        end
        
        
        % Evaluate dipole moment
        function mu = mu(obj,r)

            mu = 2.225^2 ./ ( ...
                (1 + exp(obj.b * (r{1} - obj.r0))) .* ...
                (1 + exp(obj.b * (r{2} - obj.r0))) );
            
        end
        
    end
end

