%--------------------------------------------------------------------------
%  Model potential for angular problem 
%  Dateo & Metiu, J.Chem.Phys. 95, 7392 (1991) 
%                          
%                2                     
%  V(\Theta) = \Sum  A P (cos \Theta)  
%               i=0   i i                           
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef metiu1 < pot.generic & handle
    
    properties (Access = public)
        a0
        a1
        a2
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = metiu1
        end
        
        % Initialize potential: Set/check parameters
        function init (~)
            global hamilt space
            if space.n_dim > 1
                prt.error('This potential works only for one DOF')
            end
            
            if ~isa(space.dof{1}, 'dof.legendre')
                prt.error('This potential works only for Legendre DVR')
            end
            
            if hamilt.coupling.n_eqs > 1
                prt.error('This potential works only for one Schroedinger equation')
            end
        end
    
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            prt.disp('Model potential for angular problem ')
            prt.disp('***************************************************************')
            prt.disp('                                                  ')
            prt.disp([ 'Parameter A0 : ' num2str(obj.a0) ])
            prt.disp([ 'Parameter A1 : ' num2str(obj.a1) ])
            prt.disp([ 'Parameter A2 : ' num2str(obj.a2) ])
            prt.disp (' ')
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)
            
            % Note that a call to legendre(N,x) returns all associated Legendre polynomials
            % of degree N. So we have to mask out the values for m=0 by hand and recast them
            % to the original form, which is a bit of annoying.
            leg1 = legendre(1, r{1}(:));
            leg1 = leg1(1, :)';
            leg2 = legendre(2, r{1}(:));
            leg2 = leg2(1, :)';
            
            V = obj.a0 + obj.a1 * leg1 + obj.a2 * leg2;
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            prt.error ('Code for calculation of forces still missing')
        end
        
    end
end
