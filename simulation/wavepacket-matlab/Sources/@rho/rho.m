%--------------------------------------------------------------------------
%
% Representing quantum density operators by density matricess 
% based on an eigen|energy representation
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt
%
% see the README file for license details.

classdef rho < ket & handle
    
    methods (Access = public)
        % Constructor: Setting default values
        function obj = rho ()
            
            % Inherit from constructor of "ket" superclass
            obj = obj@ket;
            
            % Setting the name of this class
            obj.string0 = 'rho';
        
        end
            
        % Logfile/console output
        function disp(~)
            prt.disp ('***************************************************************')
            prt.disp ('Matrices A, N, B and C, D and input x and output y             ')
            prt.disp ('***************************************************************')
            prt.disp ('                                                               ')
            prt.disp ('for use in a bilinear control problem                          ')
            prt.disp ('https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_abncd/')
            prt.disp ('                                                               ')
            prt.disp ('         d                                                     ')
            prt.disp ('control: -- x(t) = ( A + iu(t)N ) x(t) + iu(t)B                ')
            prt.disp ('         dt                                                    ')
            prt.disp ('                                                               ')
            prt.disp ('observe: y(t) = C x(t)                                         ')
            prt.disp ('                                                               ')
            prt.disp ('see B. Schaefer-Bung, C. Hartmann, B. Schmidt, Ch. Schuette    ')
            prt.disp ('J. Chem. Phys. 135, 014112-1-13 (2011)                         ')            
            prt.disp ('***************************************************************')
            prt.disp ('                                                       ')
            prt.disp (' coming from the quantum Liouville-von Neumann equation')
            prt.disp (' using atomic units throughout, i.e. hbar = m_e = e = 1')
            prt.disp (' for quantum systems interacting with electrical fields')
            prt.disp ('                                                       ')
            prt.disp ('  d             i                                      ')
            prt.disp (' -- rho(t) = - ---- [H - F(t) mu, rho(t)] +L [rho(t)]  ')
            prt.disp (' dt            hbar   0                     D          ')
            prt.disp ('                                                       ')
            prt.disp (' H_0 is a diagonal matrix for the unperturbed system,  ')
            prt.disp (' mu is the Hermitian matrix with dipole moments        ')
            prt.disp (' and F(t) is the electric field.                       ')
            prt.disp ('                                                       ')
            prt.disp (' with L [rho] = ... Lindblad dissipation/dephasing ... ')
            prt.disp ('       D                                               ')
            prt.disp ('                                                       ')
            prt.disp (' <C>(t) = tr( C rho(t) )                               ')
            prt.disp ('                                                       ')
        end
        
        % see separate files for the following public methods
        init_obj  ( obj, H0 )            % Initial conditions
        observe   ( obj, step )
        
    end
end

