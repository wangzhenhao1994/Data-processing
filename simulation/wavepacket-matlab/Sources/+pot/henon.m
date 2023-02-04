%--------------------------------------------------------------------------
%
% Potential energy: Henon-Heiles system (so far only in 2 dimensions)
% 
%         A    2    2         2       3                 
% V (R) = - ( R  + R ) + L ( R  R  - R /3 )             
%         2    1    2         1  2    2                 
%                                                       
% or in polar coordinates                               
%                                                       
%        A    2   1          3                          
% V(R) = - |R|  + - L |R| sin (3 theta)                 
%        2        3                                     
%                                                       
% which makes the C_3v symmetry obvious. In addition to 
% the minimum found at |R|=0, there are three saddles   
% at |R| = A/L and theta = pi/3, pi, 5*pi/3 with        
% a corresponding energy of A^3/(6*L^2).                
% Note that a classical particle can escape to infinity 
% if its energy exceeds that of the three saddle points.
% Quantum mechanical bound states can only be found     
% below the energy of the saddles. Strictly speaking,   
% these states are resonances because they can decay    
% by tunneling through the barriers.                    
%                                                       
% M.J.Davis, E.J.Heller, J. Chem. Phys. 71, 3383 (1979) 
% J.-P. Kuska / C. Herrmann, Physik Journal (Sept. 2005) 
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef henon < pot.generic & handle
    
    properties (Access = public)
        
        A           % Force constant A
        L           % Cubic parameter lambda
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = henon
        end
        
        % Initialize potential: Set/check parameters
        function init (~)
            global space
            % Check validity
            if space.n_dim ~= 2
                prt.error ('This potential only for two dimensions')
            end
            
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            prt.disp ('Potential energy: Henon-Heiles system in 2 dimensions  ')
            prt.disp ('********************************************************************')
            prt.disp ('                                                       ')
            prt.disp ( [ 'Force constant A       : ' num2str(obj.A) ] )
            prt.disp ( [ 'Cubic parameter lambda : ' num2str(obj.L) ] )
            prt.disp ( [ 'Saddle height          : ' num2str(obj.A^3/(6*obj.L^2)) ] )
            prt.disp ( [ 'Radius of triangle     : ' num2str(obj.A  /   obj.L   ) ] )
            prt.disp (' ')
        end
        
        % Evaluate potential energy functions
        function V = V(obj,r)
            V = + obj.A/2 ...
                * ( r{1}.^2 + r{2}.^2 ) ...
                + obj.L ...
                * ( r{1}.^2 .* r{2} - r{2}.^3  / 3);
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            prt.error ('Code for calculation of forces still missing')
        end
        
        
    end
end








