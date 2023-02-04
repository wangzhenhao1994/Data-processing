%--------------------------------------------------------------------------
% LEPS potential for atom/diatom A + BC reactive scattering 
%
% LEPS = London, Eyring, Polanyi, Sato             
% First two bond coordinates are assumed to        
% correspond to the first two degrees of freedom   
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%
% see the README file for license details.

classdef leps < pot.generic & handle

    properties (Access = public)
        
        d_e         % dissociation energy
        r_e         % equilibrium position
        alf         % range parameter
        s_p         % Sato parameter
        ang         % bending angle
        
    end
    
    methods (Access = public)
        
       % Constructor: Set default values
        function obj = leps
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            global space
            if space.n_dim ~= 2
                prt.error( 'This surface is for 2 dimensions.')
            end
            if length(obj.d_e) ~= 3
                prt.error( 'Wrong number of LEPS parameters')
            end
            if length(obj.r_e) ~= 3
                prt.error( 'Wrong number of LEPS parameters')
            end
            if length(obj.alf) ~= 3
                prt.error( 'Wrong number of LEPS parameters')
            end
             if length(obj.s_p) ~= 3
                prt.error( 'Wrong number of LEPS parameters')
            end
       end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            prt.disp( 'Fixed angle atom/diatom A+BC reactive scattering  ')
            prt.disp( '***************************************************************')
            prt.disp( '                                                  ')
            prt.disp( ['Bending angle       : ' num2str(obj.ang)] )
            for n = 1:3
                prt.disp( ' ' )
                prt.disp( ['Bond length #       : ' int2str(        n )] )
                prt.disp( ['dissociation energy : ' num2str(obj.d_e(n))] )
                prt.disp( ['equilibrium distance: ' num2str(obj.r_e(n))] )
                prt.disp( ['range parameter     : ' num2str(obj.alf(n))] )
                prt.disp( ['Sato parameter      : ' num2str(obj.s_p(n))] )
            end
            prt.disp (' ')
        end
        
        % Evaluate potential energy functions
        function V = V(obj,r)
            
            % Interatomic distances
            R{1} = r{1};
            R{2} = r{2};
            R{3} = sqrt(R{1}.^2 + R{2}.^2 - 2*R{1}.*R{2}*cos(obj.ang));
            
            % Coulomb and exchange integrals
            Q = cell(3,1);
            J = cell(3,1);
            for n = 1:3
                Q{n} = obj.d_e(n) / (4 + 4*obj.s_p(n))  * ( (3 + obj.s_p(n)) ...
                    * exp(-2*obj.alf(n)*(R{n} - obj.r_e(n))) ...
                    - (2 + 6*obj.s_p(n)) * exp(-obj.alf(n)*(R{n} - obj.r_e(n))) );
                
                J{n} = obj.d_e(n) / (4 + 4*obj.s_p(n)) * ( (1 + 3*obj.s_p(n)) ...
                    * exp(-2*obj.alf(n)*(R{n} - obj.r_e(n))) ...
                    - (6 + 2*obj.s_p(n)) * exp(-obj.alf(n)*(R{n} - obj.r_e(n))) );
            end
            
            % Get the complete LEPS surface
            V = Q{1} + Q{2} + Q{3} - ...
                sqrt( 0.5 * ((J{1} - J{2}).^2 + (J{1} - J{3}).^2 + (J{2} - J{3}).^2) );
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            prt.error ('Code for calculation of forces still missing')
        end
        
    end
end
