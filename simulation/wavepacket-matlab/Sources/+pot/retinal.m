% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef retinal < pot.generic & handle
    
    properties (Access = public)
        
        shift       % Vertical shift  
        barrier     % Barrier height  
        omega       % Frequency       
        kappa       % Linear term     
        
        xc          % Coupling coordinate  
        lambda      % Coupling constant 
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = retinal
        end
        
        % Initialize potential: Set/check parameters
        function init (~)
            
            global hamilt space
            
            if space.n_dim > 2
                prt.error ('This potential is only for 1 or 2 dimensions')
            end
            
            if hamilt.coupling.n_eqs ~= 2
                prt.error ('This potential is only for 2 states')
            end
            
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            global space
            disp@pot.generic (obj)
            prt.disp ('Retinal isomerization')
            prt.disp ('S.Hahn, G.Stock, JPC B 104(6),1149 (2000)')
            prt.disp ('***************************************************************')
            prt.disp ( ' ' )
            if obj.row==obj.col
                prt.disp ( [ 'Vertical shift      : ' num2str(obj.shift)   ] )
                prt.disp ( [ 'Barrier height      : ' num2str(obj.barrier) ] )
                prt.disp ( [ 'Frequency           : ' num2str(obj.omega)   ] )
                prt.disp ( [ 'Linear term         : ' num2str(obj.kappa)   ] )
            else
                prt.disp ( [ 'Coupling constant   : ' num2str(obj.lambda)  ] )
            end
            if space.n_dim == 1
                prt.disp ( [ 'Coupling coordinate : ' num2str(obj.xc)      ] )
            end
            prt.disp (' ')
        end
        
        % Evaluate potential energy functions
        function V = V(obj,r)
            global space
            if space.n_dim == 1
                coupling = obj.xc * ones(size(r{1}));
            else
                coupling = r{2};
            end
            
            if obj.row==obj.col
                V = + obj.shift ...
                    + obj.barrier * (1 - cos(r{1})) / 2 ...
                    + obj.omega * coupling.^2/2 ...
                    + obj.kappa * coupling; 
            else
                V = obj.lambda * coupling;
            end  
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            prt.error ('Code for calculation of forces still missing')
        end
        
        
    end
end


