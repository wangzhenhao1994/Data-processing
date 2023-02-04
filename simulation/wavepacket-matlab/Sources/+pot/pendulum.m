%--------------------------------------------------------------------------
% Potential energy for a generalized pendulum
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2015 Burkhard Schmidt
%
% see the README file for license details.

classdef pendulum < pot.generic & handle
    properties (Access = public)
                
        xi          % prefactor of csc^2
        chi         % prefactor of cot*csc
        eta         % prefactor of cos
        zeta        % prefactor of cos^2
        v_0         % energy shift
        
    end
        
    methods (Access = public)
    
        % Constructor: Set default values
        function obj = pendulum
            obj.xi   = 0;
            obj.chi  = 0;
            obj.eta  = 0;
            obj.zeta = 0;
            obj.v_0  = 0;
        end
        
        % Initialize potential: Set/check parameters
        function init (~)
            global space
            if space.n_dim ~= 1
                prt.error ('This potential only for one dimension')
            end
         end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            prt.disp ('Quasi-exactly solvable generalized pendulum')
            prt.disp ('***************************************************************')
            prt.disp ('                                                  ')
            prt.disp (' V = xi*csc^2+chi*cot*csc-eta*cos-zeta*cos^2+v_0  ')
            prt.disp (' see work by B. Schmidt and B. Friedrich (2014-17)')
            prt.disp ('                                                  ')
            prt.disp ('Important note on usage for planar pendulum       ')
            prt.disp ('    0 <= theta <= 2*pi (or also 4*pi !)           ')
            prt.disp ('Important note on usage for spherical pendulum    ')
            prt.disp ('    0 <= theta <= pi  and xi = m^2 - 1/4          ')
            prt.disp ('    where m is the azimuthal quantum number       ')
            prt.disp (' ')
            prt.disp ( [ 'csc^2 term  (xi)   : ' num2str(obj.xi)   ] )
            prt.disp ( [ 'cot*csc term (chi) : ' num2str(obj.chi)  ] )
            prt.disp ( [ 'Orientation (eta)  : ' num2str(obj.eta)  ] )
            prt.disp ( [ 'Alignment (zeta)   : ' num2str(obj.zeta) ] )
            prt.disp ( [ 'Energy shift (v_0) : ' num2str(obj.v_0)  ] )
            prt.disp (' ')
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)
            
            % Truncate csc^2 = 1/sin^2 to prevent singularity
            my_eps = 1e-5;
            x = r{1};
            csc2 = zeros(size(x));
            for ii=1:length(x)
                if abs(x(ii)) < my_eps || abs(x(ii)-pi) < my_eps
                    csc2(ii) = 1/sin(my_eps)^2;
                else
                    csc2(ii) = 1/sin(x(ii))^2;
                end
            end
            
            % Evaluate potential energy function
            V = ...
                + obj.xi*csc2 ...
                + obj.chi*csc2.*cos(x) ...
                - obj.eta*cos(x) ...
                - obj.zeta*cos(x).^2 ...
                + obj.v_0;
            
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            prt.error ('Code for calculation of forces still missing')
        end
        
        
        
    end
end








