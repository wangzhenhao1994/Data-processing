%--------------------------------------------------------------------------
%
% Wave functions represented on DVR / FBR grids
% for use in fully quantum-mechanical dynamics
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-.... Burkhard Schmidt
%
% see the README file for license details.

classdef wave < generic & handle
    
    properties (Access = public)
        
        dvr         % Discrete variable representation (cell vector for coupled states)
        adi         % Backup copy for adi=>dia transformation
        dia         % Backup copy for dia=>adi transformation
        fbr         % Finite basis representation (cell vector for coupled states)
        ini         % Initial wavefunction (DVR, cell vector)
        new         % New wavefunction (DVR, cell vector)
        old         % Old wavefunction (DVR, cell vector)
        sum         % Sum of wavefunctions (DVR, Chebychev only, cell vector)
        mom         % Momentum operator acting on wavefunction (DVR, single channel)
        
        bound       % Bound states; only for relaxation with 'cheby_imag'
        
        redu        % Reduced densities (for plots only)
        wig         % Wigner transform (for plots only, cell vector)
        wig_max     % Maximum of Wigner transform
        
        M_ham       % Matrix (eigen) representation: Hamiltonian 
        M_amo       % Matrix (eigen) representation: Additional mult. op.'s
        M_dip       % Matrix (eigen) representation: Dipole moment
        M_pol       % Matrix (eigen) representation: Polaritability
        M_sbc       % Matrix (eigen) representation: System-bath coupling
        M_mat       % Matrix (eigen) representation: Observables
        M_vec       % Vector (eigen) representation: Observables
        M_lab       % Labels of observables
        M_obs       % Types of observables
        
    end
    
    methods (Access = public)
        
        % Constructor: Setting default values
        function obj = wave ()
            
            % Inherit from constructor of generic superclass
            obj = obj@generic;
            
            % Setting the name of this class
            obj.string0 = 'wave';
            
            % Cell array with bound states; only for 'cheby_imag' 
            obj.bound = {};
        
        end
       
        
        % Logfile/console output
        function disp(~)
            prt.disp (' ')
            prt.disp ('***************************************************************')
            prt.disp ('Solving (coupled) Schroedinger equations by DVR/FBR techniques')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp ('For quantum systems interacting with electrical fields              ')
            prt.disp ('(using semiclassical dipole approximation)                          ')
            prt.disp ('using partial differential equations (pde) solvers                  ')
            prt.disp ('Atomic units (hbar = m_e = e = 1) are used throughout.              ')
            prt.disp ('https://sf.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_propa ')
            prt.disp ('                                                                    ')
            prt.disp ('H = T + V -iW -F*\mu -F^2/2*\alpha                                  ')
            prt.disp ('                                                                    ')
            prt.disp ('with T = T(-i d/dR)       Kinetic energy    (scalar)                ')
            prt.disp ('with V = V(R)             Potential energy  (matrix)                ')
            prt.disp ('                                                                    ')
            prt.disp ('TDSE (qm_propa) only:                                               ')
            prt.disp ('with W = W(R)             Negative imaginary potential (vector)     ')
            prt.disp ('with F = F(t)             Electrical field  (scalar, time-dependent)')
            prt.disp ('with \mu = \mu(R)         Dipole moment     (matrix)                ')
            prt.disp ('with \alpha = \alpha(R)   Polarizability    (vector)                ')
            prt.disp (' ') 
        end
        
        % How much memory needed to save a "wave" object
        function out = memory_size (~)
            global hamilt space
            out = numel(space.dvr{1}) * 16 * hamilt.coupling.n_eqs;
        end
        
        % More methods: see separate files
        init_obj  ( obj )                % Initial conditions
        propagate ( obj, step )          % Propagation
        eigen     ( obj, step )          % Eigenfunctions of Hamiltonian
        observe   ( obj, step )          % Expectation values / uncertainties
        init_ham  ( obj, e, norm )       % Application of Hamiltonian
        apply_ham ( obj, e, norm )       % Application of Hamiltonian
        adiabatic ( obj, step, direction ) % Adiabatic<=>diabatic transformation
        diabatic  ( obj )                % Adiabatic=>diabatic (initial only)
        save      ( obj, step )          % Saving dynamics to file(s)
        load      ( obj, step )          % Loading dynamics from file(s)
        save_0    ( obj       )          % Saving general info to file
        load_0    ( obj, choice )        % Loading general from file
        
    end
    
    methods (Static)
        
        wig = wigner (dvr)
        
        psi_norm = normalize ( psi_in             )
        retval   =    braket ( bra,           ket )
        retval   =  sandwich ( bra, operator, ket )
        
    end
    
end

