%--------------------------------------------------------------------------
%
% Berlin single crossing example:   Horenko, Salzmann, Schmidt, Schuette 
%
%             ( A/R   C   ) 
% V_dia (R) = (           )
%             ( C    BR^2 )
%
% The Journal of Chemical Physics 117, 11075 (2002)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef single < pot.generic & handle
    
    properties (Access = public)
        
        A           % Parameter of repulsion 
        B           % Parameter of attraction
        C           % Constant coupling      

    end
    
    methods (Access = public)
        
       % Constructor: Set default values
        function obj = single
            
            obj.A = 1.0;
            obj.B = 1.0;
            obj.C = 0.1;
            
        end
        
        % Initialize potential: Set/check parameters
        function init (~)
            
            global hamilt space
                                    
            if hamilt.coupling.n_eqs ~= 2
                prt.error ('This potential is only for 2 channels')
            end  
            
            if space.n_dim ~= 1
                prt.error ('This potential is only for 1 dimension')
            end
            
        end
                
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            if obj.row==1 && obj.col==1
                prt.disp ('Single crossing example:  J. Chemical Physics 117, 11075 (2002)')
                prt.disp ('***************************************************************')
                prt.disp (' ')
                prt.disp ( [ 'Parameter of repulsion  A : ' num2str(obj.A) ] )
                prt.disp ( [ 'Parameter of attraction B : ' num2str(obj.B) ] )
                prt.disp ( [ 'Constant coupling       C : ' num2str(obj.C) ] )
            else
                prt.disp ('Same as above')
                prt.disp ('***************************************************************')
            end
            prt.disp (' ')
        end
        
        % Evaluate potential energy functions
        function V = V(obj,r)
            if obj.row==1 && obj.col==1
                V = obj.A ./ r{1};
            elseif obj.row==2 && obj.col==2
                V = obj.B * r{1}.^2;
            elseif obj.row==1 && obj.col==2
                V = obj.C * ones(size(r{1}));
            else
                prt.error ('Invalid choice of diabatic row/column indices')
            end
        end 
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            if obj.row==1 && obj.col==1
                F{1} = obj.A ./ r{1}.^2;
            elseif obj.row==2 && obj.col==2
                F{1} = - 2 * obj.B * r{1};
            elseif obj.row==1 && obj.col==2
                F{1} = zeros(size(r{1}));
            else
                prt.error ('Invalid choice of diabatic row/column indices')
            end
        end
 
    end
    
end


