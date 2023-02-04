%--------------------------------------------------------------------------
%
% Dual crossing example
%
% John C. Tully  
% Journal of Chemical Phyics 93, 1061 (1990)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef tully2 < pot.generic & handle
    
    properties (Access = public)
        
        A
        B
        C
        D
        E
                      
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = tully2
            obj.A = 0.10;
            obj.B = 0.28;
            obj.C = 0.015;
            obj.D = 0.06;
            obj.E = 0.05;
        end
        
        % Initialize potential: Set/check parameters
        function init(~)
            
            global hamilt space
            
            if space.n_dim ~= 1
                prt.error ('This potential is only for 1 dimension')
            end
            
            if hamilt.coupling.n_eqs ~= 2
                prt.error ('This potential is only for 2 states')
            end
            
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            if obj.row==1 && obj.col==1
                prt.disp ('Dual crossing example: J. of Chemical Phyics 93, 1061 (1990)   ')
                prt.disp ('***************************************************************')
                prt.disp (' ')
                prt.disp ( [ 'Tully parameter A : ' num2str(obj.A) ] )
                prt.disp ( [ 'Tully parameter B : ' num2str(obj.B) ] )
                prt.disp ( [ 'Tully parameter C : ' num2str(obj.C) ] )
                prt.disp ( [ 'Tully parameter D : ' num2str(obj.D) ] )
                prt.disp ( [ 'Tully parameter E : ' num2str(obj.E) ] )
            else
                prt.disp ('Same as above')
                prt.disp ('***************************************************************')
            end
            prt.disp (' ')
        end
        
        % Evaluate potential energy functions
        function V = V(obj,r)
            if obj.row==1 && obj.col==1
                V = + zeros(size(r{1}));
            elseif obj.row==2 && obj.col==2
                V = - obj.A * exp ( - obj.B * r{1}.^2 ) + obj.E;
            elseif obj.row==1 && obj.col==2
                V = + obj.C * exp ( - obj.D * r{1}.^2 );
            end
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            if obj.row==1 && obj.col==1
                F{1} = + zeros(size(r{1}));
            elseif obj.row==2 && obj.col==2
                F{1} = - 2 * obj.A * obj.B * r{1} .* exp ( - obj.B * r{1}.^2 );
            elseif obj.row==1 && obj.col==2
                F{1} = + 2 * obj.C * obj.D * r{1} .* exp ( - obj.D * r{1}.^2 );
            end
        end
        
    end

end



