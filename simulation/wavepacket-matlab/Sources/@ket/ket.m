%--------------------------------------------------------------------------
%
% Representing quantum states by ket vectors 
% based on an eigen|energy representation
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt
%
% see the README file for license details.

classdef ket < generic & handle
    
    properties (Access = public)
        
        x           % input (state vector)
        y           % output (observables)
        
        x_equilib   % equilibrium state vector 
        y_equilib   % equilibrium output vector 
        
        x_initial   % equilibrium state vector 
        y_initial   % equilibrium output vector 
        
        x_label     % labeling state compontents 
        y_label     % labeling output components 
        
        A           % matrix;   from energies 
        B           % vectors;  from transition dipole moments
        N           % matrices; from transition dipole moments
        C           % vectors;  from (linear) observables
        D           % matrices; from (quadratic) observables
        Q           % whether linear observable should be squared
        S           % balancing transformation
        T           % balancing transformation
        
        norm        % Norm of state vector
        energy      % Expectation value of energy
        
        label       % for plot labels
        title       % for plot titles
        
        len_CD      % length of cell vectors C or D (qm_optimal only)
        
    end
    
    methods (Access = public)
            
        % Constructor: Setting default values
        function obj = ket ()
            
            % Inherit from constructor of generic superclass
            obj = obj@generic;
            
            % Setting the name of this class
            obj.string0 = 'ket';
        
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
            prt.disp ('                 T                                             ')
            prt.disp ('observe: y(t) = x (t) D x(t)                                   ')
            prt.disp ('                                                               ')
            prt.disp ('see B. Schaefer-Bung, C. Hartmann, B. Schmidt, Ch. Schuette    ')
            prt.disp ('J. Chem. Phys. 135, 014112-1-13 (2011)                         ')            
            prt.disp ('***************************************************************')
            prt.disp ('                                                       ')
            prt.disp (' coming from the time-dependent Schoedinger equation   ')
            prt.disp (' using atomic units throughout, i.e. hbar = m_e = e = 1')
            prt.disp (' for quantum systems interacting with electrical fields')
            prt.disp ('                                                       ')
            prt.disp ('    d                                                  ')
            prt.disp (' i -- |psi(t)> = ( H  - F(t) mu ) |psi(t)>             ')
            prt.disp ('   dt               0                                  ')
            prt.disp ('                                                       ')
            prt.disp (' H_0 is a diagonal matrix for the unperturbed system,  ')
            prt.disp (' mu is the Hermitian matrix with dipole moments        ')
            prt.disp (' and F(t) is the electric field.                       ')
            prt.disp ('                                                       ')
            prt.disp ('  <D>(t)= <psi(t)| D |psi(t)>                          ')
            prt.disp ('                                                       ')
            
        end
        
        % How much memory needed to save a "ket" object
        function out = memory_size (obj)
            out = length(obj.x) * 16;
        end

        % see separate files for the following public methods
        init_obj     ( obj, H0 )         % Initial conditions
        propagate    ( obj, step )       % Solving the TDSE using ODE methods
        observe      ( obj, step )       % Expectation values from CD
        init_ham     ( obj )             % Initialization of ABNCD
        apply_ham    ( obj )             % Application of ABN
        save         ( obj, step )       % Saving dynamics to file(s)
        load         ( obj, step )       % Loading dynamics from file(s)
        save_0       ( obj       )       % Saving general info to file
        load_0       ( obj, choice )     % Loading general from file

        spectrum_A   ( obj,choice )      % plot spectrum of A matrix
        check_stable ( obj, label )      % check stability of A evolution
        adiabatic    ( ~, ~, ~ )         % Adiabatic<=>diabatic transformation


    end
end

