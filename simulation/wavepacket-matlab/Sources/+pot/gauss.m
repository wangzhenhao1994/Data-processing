%--------------------------------------------------------------------------
%
% Potential energy function:
% Gaussian bell-shaped function 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2020-20xy Burkhard Schmidt
%
% see the README file for license details.

classdef gauss < pot.generic & handle
    
    properties (Access = public)
        height      % Height of Gaussian
        pos_0       % Central position of Gaussian
        width       % Width of Gaussian
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = gauss
            obj.height = 1;
            obj.pos_0  = 0;
            obj.width  = 1;
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            
            global space
            
            if ~isscalar(obj.height)
                prt.error ('Incompatible dimensionality for height')
            end
            
            if length(obj.pos_0)~=space.n_dim
                prt.error ('Incompatible dimensionality for pos_0')
            end
            
            if length(obj.width)~=space.n_dim
                prt.error ('Incompatible dimensionality for width')
            end
                        
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic(obj)
            prt.disp ( 'Gaussian bell-shaped potential energy function in N-dim' )
            prt.disp ( '***************************************************************' )
            prt.disp ( '                                                      ' )
            prt.disp ( '               N      [   ( Ri-R0i )^2 ]              ' )
            prt.disp ( ' V (R) = H * Prod exp [ - (--------)   ]              ' )
            prt.disp ( '              i=1     [   (  2*Wi  )   ]              ' )
            prt.disp ( '                                                      ' )
            prt.disp ( 'where the product extends over all spatial dimensions ' )
            prt.disp ( '***************************************************************' )
            prt.disp (   ' ')
            prt.disp ( [ 'Height of Gaussian         H : ' num2str(obj.height) ] )
            prt.disp ( [ 'Mean value position       R0 : ' num2str(obj.pos_0 ) ] )
            prt.disp ( [ 'Position uncertainty       W : ' num2str(obj.width ) ] )
            prt.disp (' ')
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)
            global space
            
            V = obj.height * ones ( size(r{1}) );
            
            % Tensor product of one-dimensional Gaussians
            for k = 1:space.n_dim
                V = V .* exp (  -((r{k}-obj.pos_0(k)) / (2*obj.width(k))).^2  );
            end
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            
            F = cell(size(r));
            
            V = obj.height * ones ( size(r{1}) );
            for k = 1:length(r)
                V = V .* exp (  -((r{k}-obj.pos_0(k)) / (2*obj.width(k))).^2  );
            end      
            
            for d=1:length(r)
                F{d} = V * ( r{d}-obj.pos_0(d) ) / obj.width(d);
            end
        end

    end
end


