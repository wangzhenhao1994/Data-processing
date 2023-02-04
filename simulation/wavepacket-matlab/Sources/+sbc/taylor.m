%--------------------------------------------------------------------------
% Represent the system-bath coupling as a Taylor series
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2016-2017 Burkhard Schmidt
%
% see the README file for license details.

classdef taylor < sbc.generic & handle
    
    properties (Access = public)
        
        hshift      % horizontal shift (scalar)
        vshift      % vertical shift (row vector)
        coeffs      % coefficients, i.e. derivatives
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = taylor
            obj.hshift  = [];
            obj.vshift  = 0;
            obj.coeffs  = [];
            obj.row     = 1;
            obj.col     = 1;
            obj.dvr     = [];
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            global space
            if isempty(obj.hshift)
                obj.hshift = zeros(1,space.n_dim);
            end
        end
        
        % Display potential, overloading default disp method
        function disp (obj)
            disp@sbc.generic(obj)
            prt.disp ('Taylor series (diag. in N dimensions)')
            prt.disp ('***************************************************************')
            prt.disp (' ')
        end
        
        % Evaluate system-bath coupling
        function chi = chi(obj,r)

            chi = math.taylor (...
                r, ...
                obj.hshift, ...
                obj.vshift, ...
                obj.coeffs, 1 );
        end
    end
end
