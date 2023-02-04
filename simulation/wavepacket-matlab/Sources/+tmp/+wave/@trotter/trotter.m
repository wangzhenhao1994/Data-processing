%-------------------------------------------------------------------------
%
% Operator splitting
% ==================
%
% Propagate objects of class wave (wavefunctions on grids) 
% by one substep of size time.steps.s_delta
%
% Hamiltonian
% -----------
%
%                                          F^2(t)          (       d)
%   H = V (R) - F(t) * ( D (R) + D (R) ) - ------ P(R) + T (R, -i --)
%                         p       t          2             (      dR)
%
%
% with kinetic operator T and potential V (R) energy (may be matrix-valued) 
% and where the interaction between an (optional) external electric 
% field F(t) and the dipole moment D(R) and/or the polarizability P(R)
% of the system is treated semiclassically. Here, the dipole matrix D 
% has been split into its diagonal (D_p) and offdiagonal (D_t) parts
% containing permanent and transition dipoles, respectively. In most
% cases, only one of the two is expected to play a role, depending on 
% the light frequencies under consideration. Moreover, note that F(t) 
% as well D_p(R), D_t(R), and P(R) can have two cartesian components 
% (along x,y) corresponding to different polarization directions.
%
% Lie-Trotter splitting
%
% psi(t+tau) = exp( -i*         V *tau ) 
%            * exp( +i*F(t)    *Dp*tau ) 
%            * exp( +i*F(t)    *Dt*tau ) 
%            * exp( +i*F^2(t)/2*P *tau ) 
%            * exp( -i*         T *tau )
%            * psi(t)                    + O(tau^2)
%
%
% Implementation for a single Schroedinger equation
% -------------------------------------------------
%
% The above Trotter and Strang formulae are straightforward to evaluate for
% FFT grids (plane wave DVR):
%   Operator V, and hence exp(V), are diagonal in position representation
%   Operator D, and hence exp(D), are diagonal in position representation
%   Operator T, and hence exp(T), are diagonal in momentum representation
% Hence, two FFTs have to be performed for every timestep to transform
% from position to momentum representation and back. Similar arguments
% apply for other DVR methods.
%
% Implementation for coupled Schroedinger equations
% -------------------------------------------------
%
% In case of coupled Schroeinger equations, the potential energy V as well
% as the dipole moment operator D may be represented by (real, symmetric) 
% matrices. Hence, for each spatial discretization point, the exponential 
% of these matrices has to be calculated. This can be achieved either by
% Matlab's expm function or by the method of eigenvectors and eigenvalues, 
% where the accuracy is determined by the condition of the matrix.
% 
% For the case of two equations, the matrix exponential can be calculated 
% analytically, see  Eq. (11.204) on page 316 of David Tannor's book. 
% 
%     ( alpha   beta )
% V = (              )
%     ( beta^* gamma )
%
% delta = (alfa-gamma)/2
%   eta = (alfa+gamma)/2
%   rho = sqrt(delta^2+beta^2)
%
% expm(-iV*tau) 
%   = exp(-i*eta*tau) *
%   * [              cos(tau*rho)*I 
%       -i*delta/rho*sin(tau*rho)*S3
%       -i* beta/rho*sin(tau*rho)*S1  ]
%
% where I=(1 0; 0 1), S1 = (0 1; 1 0), S3 = (1 0; 0 -1)
%
% 
% References
% ==========
%
% Original version: 
%     J. A. Fleck et al.,                  Appl. Phys.,   10,  129 (1976)
%     M. D. Feit, J. A. Fleck, A. Steiger, J. Comp. Phys. 47,  412 (1982) 
% Coupled equations: 
%     J. Alvarellos and H. Metiu,          J. Chem. Phys. 88, 4957 (1988)
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

classdef trotter < handle
    
    properties (Access = public)
        frac        % fraction: One (Trotter) or one half (Strang)
    end
    
    % For future versions: add these private properties (from global hamilt)
        % pot_couple
        % pot_expo{m,n}
        % dip_prod{p}{m}
        % pol_prod{p,q}{m}
        % dip_trans{p}
        % dip_eig_vals{p}{m}
        % dip_eig_vecs{p}{m,n}
        % dip_diagonal{p}{m}
    % end
    
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = trotter
            obj.frac = 1;
        end
              
        % Display propagator, overloading default disp method
        function disp(~)

            prt.disp ('Lie-Trotter splitting (first order)')
            prt.disp ('***************************************************************')
            
        end

        % Initialize propagator
        function init (obj,~)
            global space time hamilt
                        
            % Propagator for kinetic energy (all but the last one, see below)
            for k = 1:space.n_dim-1
                init_kin(space.dof{k}, obj.frac, false);
            end
            
            % External kinetic energies; Since they are typically the most
            % time-consuming operators, we put them in the middle of the split operator
            % If none are specified, take the last grid-internal kinetic energy
            % in the middle.
            if isfield(hamilt, 'kin')
                for n = 1:length(hamilt.kin)-1
                    init_kin(hamilt.kin{n}, obj.frac, false);
                end
                init_kin(hamilt.kin{end}, 1, false);
                init_kin(space.dof{end}, obj.frac, false);
            else
                init_kin(space.dof{end}, 1, false);
            end
            
            % Propagator for the potential energy
            pot_init (obj)
            
            % Propagator for (permanent and/or transition) dipoles, polarizabilities
            if isfield(time,'pulse')
                perm_init  (obj);
                trans_init (obj);
                pol_init   (obj);
            end
            
        end
        
        % Perform propagation
        function propa (obj, psi, k)
            global space time hamilt
            
            % Potential energy: Full (Trotter) or half (Strang) sub-step
            pot_propa (obj, psi);
            
            % Permanent/transition dipoles, polarizabilities: 
            % Full (Trotter) or half (Strang) sub-step
            if isfield(time,'pulse')
                e = [ ...
                    time.efield.grid{1}(k+time.steps.offset) ...
                    time.efield.grid{2}(k+time.steps.offset)];
                perm_propa  (obj, psi, e);
                trans_propa (obj, psi, e, k==1);
                pol_propa   (obj, psi, e);
            end
            
            % Kinetic energy: Propagate each DOF
            for l = 1:space.n_dim
                kinetic_exp(space.dof{l}, psi);
            end
            
            % Propagate each external kinetic energy
            if isfield(hamilt, 'kin')
                for n = 1:length(hamilt.kin)
                    kinetic_exp(hamilt.kin{n}, psi);
                end
            end
            
        end
        
    end
    
    methods (Access = public) % because these methods are inherited by strang
        
        pot_init   (obj)
        pol_init   (obj)
        trans_init (obj)
        perm_init  (obj)
        
        pot_propa   (obj, psi)
        pol_propa   (obj, psi, e);
        trans_propa (obj, psi, e, recalc);
        perm_propa  (obj, psi, e);
        
    end
end
        
    
