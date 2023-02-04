%------------------------------------------------------------------------------
%
% This function creates associated Legendre polynomials in cos Theta
%
% Apart from normalization and missing azimuthal functions,
% these functions are identical to spherical harmonics.
%
% Optionally, dividing wavefunction by sqrt(sin(Theta))
% see Eq. (4) in doi:10.1103/PhysRev.A.91.022111
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 Burkhard Schmidt
%
% see the README file for license details.

classdef spherical < init.generic & handle
    
    properties (Access = public)
        l           % Angular momentum quantum number
        m           % Azimuthal quantum number
        sqst        % Toggle division by sqrt(sin(Theta))

    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = spherical
            obj.l = 0;
            obj.m = 0;
            obj.sqst = false;
        end
        
        % Initialize spherical functions: Set/check parameters
        function init (obj)
            
            if obj.l<0
                prt.error ('quantum number l must not be negative')
            end
            if obj.m<0
                prt.error ('quantum number m must not be negative')
            end
            if obj.m>obj.l
                prt.error ('quantum number m must not exceed l')
            end
            
        end

        % Display pendular function, overloading default disp method
        function disp(obj)
            
            prt.disp ('Spherical harmonic wavefunction')
            prt.disp ('***************************************************************')
            prt.disp ('   ' )
            prt.disp ( ['Angular momentum quantum number l : ' int2str(obj.l)] )
            prt.disp ( ['Azimuthal quantum number        m : ' int2str(obj.m)] )
            prt.disp ( ['Divide by sqrt(sin(Theta))        : ' int2str(obj.sqst)] )
            
        end

        % Evaluate wave function on a grid  ==> Q/M propagation
        function wave (obj)
            global space
            
            X = cos ( space.dvr{obj.dof} );
            Y = legendre (obj.l,X);
            obj.dvr = Y (obj.m+1,:)';
            
            if obj.sqst
                obj.dvr = obj.dvr .* sqrt( sin ( space.dvr{obj.dof} ) );
            end
        end
        
        % Sample phase space density ==> Q/C propagation
        function traj (~, ~)
            prt.error ('Code for phase space sampling still missing')
        end
        
    end
end
