%--------------------------------------------------------------------------
% Ultrafast S2-S1 internal conversion of pyrazine  
% Two state / three mode model exhibits a conical  
% intersection along a straight line located in    
% the (R1,R2) plane. Note that this intersection   
% is not(!) due to a Jahn-Teller effect but is     
% only partially symmetry induced!                 
%                                                  
% Two non-degenerate electronic states (D_2h group)
%     S1 (B3u)                                     
%     S2 (B2u)                                     
% Two tuning modes (R1,R2), one coupling mode (R3) 
%     R1 = mode 01  (A_g, totally symmetric)        
%     R2 = mode 06a (A_g, totally symmetric)       
%     R3 = mode 10a (B1g symmetry)                  
%                                                  
%              3  omega_k   ( R_k^2   0  )         
%  V    (R) = Sum ------- * (            )         
%   dia       k=1   2       (  0   R_k^2 )         
%                                                  
%              2  ( kappa1_i R_i      0        )   
%           + Sum (                            )   
%             i=1 (     0         kappa2_i R_i )   
%                                                  
%              3             (  0   R_j )          
%           + Sum lambda12 * (          )          
%             j=3            ( R_j   0  )          
%                                                  
%             ( E1   0  )                          
%           + (         )                          
%             ( 0    E2 )                          
%                                                  
%  Note that the summations over kappa and lambda  
%  extend over tuning and coupling modes, respect.,
%  while the omega summation extends over all modes
%                                                  
%  Schneider, Domcke,         CPL 150, 235 (1988)  
%  Schneider, Domcke, Köppel, JCP 92, 1045 (1990)  
%
%--------------------------------------------------------------------------


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef pyrazine < pot.generic & handle
    
    properties (Access = public)
        
        omega       % Force constant omega
        kappa       % Linear parameters kappa
        lambda      % Linear coupling lambda
        energ       % Vertical excitation E 
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = pyrazine
            obj.omega  = zeros(3,1);
            obj.kappa  = zeros(2,1);
            obj.lambda = 0;
            obj.energ  = 0;
        end
        
        % Initialize potential: Set/check parameters
        function init (~)
            
            global hamilt space
            
            if space.n_dim ~= 3
                prt.error ('This potential is only for 3 dimension')
            end
            
            if hamilt.coupling.n_eqs ~= 2
                prt.error ('This potential is only for 2 equations')
            end
            
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            if obj.row==1 && obj.col==1
                prt.disp ('Ultrafast S2-S1 internal conversion of pyrazine  ')
                prt.disp ('*************************************************')
                prt.disp ('                                                 ')
                prt.disp ( [ 'Force constants omega   : ' num2str(obj.omega) ] )
                prt.disp ( [ 'Linear parameters kappa : ' num2str(obj.kappa) ] )
                prt.disp ( [ 'Vertical excitation E   : ' num2str(obj.energ) ] )
                prt.disp ( [ 'Linear coupling lambda  : ' num2str(obj.lambda) ] )
            else
                prt.disp ('Same as above')
                prt.disp ('***************************************************************')
            end
            prt.disp (' ')
        end
        
        % Evaluate grid representation of potential energy functions
        function V = V(obj,r)
            
            if (obj.row==obj.col)
                V = + obj.omega(1)/2 * r{1}.^2 ...
                    + obj.omega(2)/2 * r{2}.^2 ...
                    + obj.omega(3)/2 * r{3}.^2 ...
                    + obj.kappa(1)   * r{1} ...
                    + obj.kappa(2)   * r{2} ...
                    + obj.energ;
            else
                V = obj.lambda * r{3};
            end
            
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            prt.error ('Code for calculation of forces still missing')
        end
 
        
    end
end
