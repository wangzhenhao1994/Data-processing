%--------------------------------------------------------------------------
%
% Additional multiplicative operator:
% Educts=Reactant versus Products in a chemical exchange reaction
%
%     A + BC -> AB + C       
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%               2004-2017 Burkhard Schmidt
%
% see the README file for license details.


classdef reaction < amo.generic & handle
    
    properties (Access = public)
        reac        % Index of  reactant channel
        prod        % Index of  product channel
        side        % Choose reactant or product side        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = reaction
            obj.reac = 1;
            obj.prod = 2;
            obj.side = 'p';
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            global space
            
            if space.n_dim < 2
                prt.error('At least two degrees of freedom required!');
            end

            if obj.reac > space.n_dim
                prt.error('Wrong choice for reactant channel');
            end

            if obj.prod > space.n_dim
                prt.error('Wrong choice for product channel');
            end

           if obj.prod == obj.reac
                prt.error('Wrong choice for reactant/product channel');
            end
            
            if isempty (obj.label)
                switch lower(obj.side)
                    case 'r'
                        obj.label = 'Reactant';
                    case 'p'
                        obj.label = 'Product';
                    otherwise
                        prt.error('Reactant/product side must be "r" or "p"')
                end
            end
            
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@amo.generic(obj)
            prt.disp ('Reactants/products for A + BC -> AB + C  reaction ')
            prt.disp ('***************************************************************')
            prt.disp ('  ')
            prt.disp(['Index of reactant distance coordinate AB :' num2str(obj.reac)])
            prt.disp(['Index of product  distance coordinate BC :' num2str(obj.prod)])
            switch lower(obj.side)
                case 'r'
                    prt.disp('Choosing reactant side')
                case 'p'
                    prt.disp('Choosing product side')
            end
            prt.disp ( [ 'Label                            : '         obj.label] )
            prt.disp (' ')
        end

        % Evaluate AMO function
        function A = A (obj,r)
           
            switch lower(obj.side)
                case 'r'
                    A = ...
                        (r{obj.reac} > r{obj.prod});
                case 'p'
                    A = ...
                        (r{obj.reac} < r{obj.prod});
            end
            
        end
    end
end
