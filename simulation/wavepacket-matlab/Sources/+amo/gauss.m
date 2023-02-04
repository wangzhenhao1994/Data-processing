%--------------------------------------------------------------------------
%
% Additional multiplicative operator:
% Gaussian bell-shaped function as projector 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef gauss < amo.generic & handle
    
    properties (Access = public)
        pos_0       % Central position of Gaussian
        width       % Width of Gaussian
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = gauss
        end
        
        % Initialize AMO: Set/check parameters
        function init (obj)
            
            global space
            
            if length(obj.pos_0)~=space.n_dim
                prt.error ('Incompatible dimensionality for pos_0 of projection')
            end
            
            if length(obj.width)~=space.n_dim
                prt.error ('Incompatible dimensionality for width of projection')
            end
            
            if isempty(obj.label)
                obj.label = 'Gaussian';
            end
            
        end
        
        % Display AMO, overloading default disp method
        function disp(obj)
            disp@amo.generic(obj)
            prt.disp ( 'Gaussian bell-shaped function in N-dim' )
            prt.disp ( '***************************************************************' )
            prt.disp ( '                                                  ' )
            prt.disp ( '           N      [   ( Ri-R0i )^2 ]              ' )
            prt.disp ( ' G (R) = Prod exp [ - (--------)   ]              ' )
            prt.disp ( '          i=1     [   (  2*Wi  )   ]              ' )
            prt.disp ( '                                                  ' )
            prt.disp ( 'where the product extends over all dimensions     ' )
            prt.disp ( '                                                  ' )
            prt.disp ( 'For example, see Eq. (49) in JCP 109, 385 (1998)  ' )
            prt.disp ( 'DOI:10.1063/1.476575 by W. Zhu and H. Rabitz      ' )
            prt.disp ( '***************************************************************' )
            prt.disp (   ' ')
            prt.disp ( [ 'Mean value position       R0 : ' num2str(obj.pos_0) ] )
            prt.disp ( [ 'Position uncertainty      W  : ' num2str(obj.width) ] )
            prt.disp ( [ 'Label                        : '         obj.label  ] )
            prt.disp (' ')
        end
        
        % Evaluate AMO function
        function A = A(obj,r)
            global space
            
            A = ones ( size(r{1}) );
            
            % Tensor product of one-dimensional Gaussians
            for k = 1:space.n_dim
                A = A / ...
                    (2*obj.width(k)*sqrt(pi)) .* ...
                    exp (  -((r{k}-obj.pos_0(k)) / (obj.width(k)*2)).^2  );
            end
        end
    end
end


