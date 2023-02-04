%--------------------------------------------------------------------------
% Trigonometric (cosine-shaped) model for dipole moment
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 Burkhard Schmidt
%
% see the README file for license details.


classdef cosine < dip.generic & handle
    
    properties (Access = public)
        
        pre         % prefactor
        exp         % exponent
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = cosine
            obj.pre = 1;
            obj.exp = 1;
        end
        
        % Initialize dipole moment: Set/check parameters
        function init (~)
            global hamilt space
            
            if hamilt.coupling.n_eqs ~= 1
                prt.error ('This dipole function is only for 1 state')
            end
            
            if space.n_dim > 1
                prt.error ('This dipole function is only for 1 dof')
            end
        end
        
        % Display dipole moment, overloading default disp method
        function disp (obj)
            disp@dip.generic(obj)
            prt.disp (' mu(Theta) = f * cos^n  Theta                     ')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp (['Prefactor f : ' num2str(obj.pre)])
            prt.disp (['Exponent  n : ' num2str(obj.exp)])
            prt.disp (' ')
        end
        
        % Evaluate dipole moment
        function mu = mu(obj,r)
            global hamilt space
            if isa (space.dof{hamilt.amo{1}.dof}, 'dof.legendre')
                mu = obj.pre *     r{1} .^obj.exp;
            else
                mu = obj.pre * cos(r{1}).^obj.exp;
            end
            
        end
    end
end
