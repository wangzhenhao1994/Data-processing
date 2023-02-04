%--------------------------------------------------------------------------
% Represent the potential energy function as a Taylor series
% Includes free particle (V=0) as a special case
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2016-2017 Burkhard Schmidt
%
% see the README file for license details.

classdef taylor < pot.generic & handle
    
    properties (Access = public)
        
        hshift      % Horizontal shift (row vector)
        vshift      % Vertical shift (scalar)
        coeffs      % Coefficients, i.e. derivatives
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = taylor
            obj.hshift  = [];
            obj.vshift  = 0;
            obj.coeffs  = [];
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            global space
            if isempty(obj.hshift)
                obj.hshift = zeros(1,space.n_dim);
            end
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            if obj.vshift==0 && isempty (obj.coeffs)
                prt.disp ('Not available ("free particle")                   ')
                prt.disp ('***************************************************************')
            else
                prt.disp ('Taylor series (diag. in N dimensions)')
                prt.disp ('***************************************************************')
                prt.disp (' ')
                if any(obj.hshift)
                    prt.disp ( [ 'Horizontal shift    : ' num2str(obj.hshift) ] )
                end
                if obj.vshift
                    prt.disp ( [ 'Vertical shift      : ' num2str(obj.vshift) ] )
                end
                if ~isempty (obj.coeffs)
                    prt.disp ( [ 'Taylor coefficients : ' num2str(obj.vshift) ] )
                end
            end
            prt.disp (' ')
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)  
            V = math.taylor (...
                r, ...
                obj.hshift, ...
                obj.vshift, ...
                obj.coeffs, 0 );
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            F = math.taylor_d (...
                r, ...
                obj.hshift, ...
                obj.coeffs );
        end
        
    end
end
