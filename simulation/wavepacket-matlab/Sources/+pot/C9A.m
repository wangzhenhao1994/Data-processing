%------------------------------------------------------------------------------
%
% Potential energy function for the C9A molecule. Since the potentials
% that are fitted from spectroscopic data are completely nuts 
% (the coeffcients by Monte et al. give _something_, but not the
% correct potentials), we use a fourth order potential in the torsion
% angle. The diabatic coupling between the excited state S1 and Sx
% potentials is given by a Gaussian.
%
% S0, S1 and "dark" Sx state. 
% Potentials taken from
% Manz et al. Z.Phys.D 34:111     
% Monte et al. JChemPhys 98:2580  
% 
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 - .... Burkhard Schmidt
%               2009 - 2011 Ulf Lorenz
%
% see the README file for license details.

classdef C9A < pot.generic & handle

    properties (Access = public)
        
        state       % electronic state: S0 or S1 / Sx
        
        V0          % Coefficient of S0, Sx
        V2          % Coefficient of S0, Sx
        V4          % Coefficient of S0
        
        coeffs      % Coefficients of S1 state
        vshift      % Shift of S1 state 
        
        C           % Coefficients of S1-Sx coupling
        a           % Coefficients of S1-Sx coupling
        phi0        % Coefficients of S1-Sx coupling
        
    end
    
    methods (Access = public)
        
       % Constructor: Set default values
        function obj = C9A (choice)
            switch choice
                case 'S0'
                case 'S1'
                case 'Sx'
                case 'S1-Sx'
                otherwise
                    prt.error ('Wrong specification of electronic states')
            end
            obj.state = choice;
        end
        
        
        % Initialize potential: Set/check parameters
        function init (~)
            
            global space
            
            if space.n_dim ~= 1
                prt.error ('This potential is for one degrees of freedom only')
            end
            
        end

        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            prt.disp ('Potential for 9-(N-carbazolyl) anthracene (C9A)')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp (['Electronic state(s) : ' obj.state])
            prt.disp (' ')
        end
        
        % Evaluate grid representation of potential energy functions
        function V = V(obj,r)
            
            switch obj.state
                case 'S0'
                    V = ...
                        + obj.V0 ...
                        + obj.V2/ 2 * (r{1}-pi/2).^2 ...
                        + obj.V4/24 * (r{1}-pi/2).^4;
                    
                case 'S1'
                        V = zeros(size(r{1}));
                        for ii = 1:length(obj.coeffs)
                            V = V ...
                                + obj.coeffs(ii)/2 ...
                                * (1 - cos(2*ii*r{1}));
                        end
                        V = - V + obj.vshift;
                        
                case 'Sx'
                        V = ...
                            + obj.V0 ...
                            + obj.V2/2 * (r{1}-pi/2).^2;
                        
                case 'S1-Sx'                        
                        V = obj.C * ( ...
                            + exp( - (r{1}-pi/2+obj.phi0).^2/(2*obj.a^2) )  ...
                            + exp( - (r{1}-pi/2-obj.phi0).^2/(2*obj.a^2) ) );
            end
            
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            prt.error ('Code for calculation of forces still missing')
        end
        
   
    end
end
