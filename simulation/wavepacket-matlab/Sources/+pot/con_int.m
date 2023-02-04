% *************************************************************************
% Generic E x e conical intersection example :     
%                                                  
% Two degenerate electronic states coupled by      
% two degenerate normal modes of vibration with    
% linear and quadratic Jahn-Teller coupling        
%                                                  
%             omega   ( X^2+Y^2   0  )             
%  V    (R) = ----- * (              )             
%   dia         2     (  0   X^2+Y^2 )             
%                                                  
%                     ( X   Y )                    
%           + kappa * (       )                    
%                     ( Y  -X )                    
%                                                  
%             gamma   ( +(X^2-Y^2)  -2*X*Y   )     
%           + ----- * (                      )     
%               2     (  -2*X*Y   -(X^2-Y^2) )     
%                                                  
%  The eigenvalues give the 2 adiabatic potentials 
%                                                  
%                    2       (      2  2           
%  E    (R) = omega R /2 +/- ( kappa  R  + ...     
%   adi                      (                     
%                                                  
%                 3                   2  4   ) 1/2 
%  + gamma*kappa*R *cos(3*phi) + gamma  R /4 )     
%                                            )     
%                                                  
%  with polar coordinates                          
%  R = ( X^2 + Y^2 )^1/2 and phi = atan(Y/X)       
%                                                  
% see e. g. 
%
% A. Viel and W. Eisfeld 
% Effects of higher order Jahn-Teller coupling on the nuclear dynamics
% J. Chem. Phys. 120, 4603  (2004) 
% https://doi.org/10.1063/1.1646371
%
% W. Eisfeld and A. Viel   
% Higher order (A+E) x e pseudo-Jahn–Teller coupling
% J. Chem. Phys. 122, 204317  (2005) 
% https://doi.org/10.1063/1.1904594
% where the (pseudo) Jahn-Teller is discussed to 6th order !
%                                                  
% *************************************************************************

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef con_int < pot.generic & handle
    
    properties (Access = public)
        
        omega       % Force constant omega
        kappa       % Linear JT coupling
        gamma       % Quadratic JT coupling
        delta       % Energy gap
        
        x_ci        % x-shift of conical intersection
        y_ci        % y-shift of conical intersection
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = con_int
            obj.omega = 1;
            obj.kappa = 1;
            obj.gamma = 0;
            obj.delta = 0;
            
            obj.x_ci  = 0;
            obj.y_ci  = 0;
        end
        
        % Initialize potential: Set/check parameters
        function init (~)
            
            global hamilt space
            
            if space.n_dim > 2
                prt.error ('This potential is only for 1 or 2 dimension')
            end
            
            if hamilt.coupling.n_eqs > 2
                prt.error ('This potential is only for 1 or 2 (coupled) channels')
            end
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            
            global space
            
            disp@pot.generic (obj)
            if obj.row==1 && obj.col==1
                
                prt.disp ('Generic E x e conical intersection example :     ')
                prt.disp ('***************************************************************')
                prt.disp ('                                                 ')
                prt.disp ( [ 'Force constant            omega : ' num2str(obj.omega) ] )
                prt.disp ( [ 'Linear JT-coupling        kappa : ' num2str(obj.kappa) ] )
                prt.disp ( [ 'Quadr. JT-coupling        gamma : ' num2str(obj.gamma) ] )
                if space.n_dim == 1
                    prt.disp ( [ 'Energy gap                delta : ' num2str(obj.delta) ] )
                end
                if obj.x_ci ~= 0
                    prt.disp ( [ 'x-shift of conical intersection : ' num2str(obj.x_ci) ] )
                end
                if obj.y_ci ~= 0
                    prt.disp ( [ 'y-shift of conical intersection : ' num2str(obj.y_ci) ] )
                end
                
            else
                prt.disp ('Same as above')
                prt.disp ('***************************************************************')
            end
            prt.disp (' ')
            
        end
        
        % Evaluate potential energy functions
        function V = V(obj,r)
            
            global hamilt space
            
            % Shifted coordinates
            x = r{1}-obj.x_ci;
            x2 = x.^2;
            if space.n_dim==2
                y = r{2}-obj.y_ci;
                y2 = y.^2;
                r2 = x2 + y2;
            elseif space.n_dim==1
                r2 = x2;
            end
            
            eps1 = obj.omega/2 * r2 + obj.kappa * x;
            eps2 = obj.omega/2 * r2 - obj.kappa * x;
            if space.n_dim==2
                beta  =     + obj.kappa   *     y;
                beta = beta - obj.gamma   * (x.*y);
                eps1 = eps1 + obj.gamma/2 * (x2-y2);
                eps2 = eps2 - obj.gamma/2 * (x2-y2);
            elseif space.n_dim==1
                beta  =     + obj.gamma;
                eps1 = eps1 + obj.delta/2;
                eps2 = eps2 - obj.delta/2;
            end
            
            % Coupled two state problem
            if hamilt.coupling.n_eqs==2
                if obj.row==1 && obj.col==1
                    V = eps1;
                end
                if obj.row==2 && obj.col==2
                    V = eps2;
                end
                if obj.row==1 && obj.col==2
                    V = beta;
                end
                
            % Only lower adiabatic state: Mexican hat
            elseif hamilt.coupling.n_eqs==1
                eta = (eps1+eps2)/2;
                dlt = (eps1-eps2)/2;
                rho = sqrt(dlt.^2+beta.^2);
                V = eta - rho;
            end
            
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            global hamilt space
            if obj.gamma~=0
                prt.error ('Code for quadratic coupling missing')
            end
            if space.n_dim==1
                prt.error ('Code for one dimension missing')
            end
            
            % Shifted coordinates
            x = r{1}-obj.x_ci;
            y = r{2}-obj.y_ci;
            
            % Coupled two state problem
            if hamilt.coupling.n_eqs==2
                
                if obj.row==1 && obj.col==1
                    F{1} = - obj.omega * x - obj.kappa;
                    F{2} = - obj.omega * y;
                end
                if obj.row==2 && obj.col==2
                    F{1} = - obj.omega * x + obj.kappa;
                    F{2} = - obj.omega * y;
                end
                if obj.row==1 && obj.col==2
                    F{1} =              zeros(size(x));
                    F{2} = - obj.kappa * ones(size(x));
                end
                
            % Only lower adiabatic state: Mexican hat
            elseif hamilt.coupling.n_eqs==1
                rr = sqrt(x.^2 + y.^2);
                F{1} = - obj.omega * x + obj.kappa * x ./ rr;
                F{2} = - obj.omega * y + obj.kappa * y ./ rr;
            end
            
        end        
        
    end
end
