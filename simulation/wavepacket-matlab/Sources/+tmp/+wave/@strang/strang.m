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
% Strang-Marchuk splitting
%
% psi(t+tau) = exp( -i*         V *tau/2 ) 
%            * exp( +i*F(t)    *Dp*tau/2 ) 
%            * exp( +i*F(t)    *Dt*tau/2 ) 
%            * exp( +i*F^2(t)/2*P *tau/2 ) 
%            * exp( -i*         T *tau   ) 
%            * exp( +i*F^2(t)/2*P *tau/2 ) 
%            * exp( +i*F(t)    *Dt*tau/2 )
%            * exp( +i*F(t)    *Dp*tau/2 )
%            * exp( -i*         V *tau/2 ) 
%            * psi(t)                    + O(tau^3)
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

classdef strang < tmp.wave.trotter & handle
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = strang
            obj.frac = 1/2;
        end
              
        % Display propagator, overloading default disp method
        function disp(~)

            prt.disp ('Strang-Marchuk splitting (second order)')
            prt.disp ('***************************************************************')
                       
        end

        
        % Perform propagation
        function propa (obj, psi, k)
            global space time hamilt
            
            % Same as in Trotter, but with frac=1/2
            propa@tmp.wave.trotter (obj, psi, k)

            % Inverse order of the external kinetic energies
            % And apply the last grid kinetic operator as well.
            if isfield(hamilt, 'kin')
                for n = length(hamilt.kin)-1:-1:1
                    kinetic_exp(hamilt.kin{n}, psi);
                end
                kinetic_exp(space.dof{end}, psi);
            end
            
            % Apply all but the last grid kinetic energy operator
            % (it has either been applied before, or it was the
            % innermost split operator).
            for l = space.n_dim-1:-1:1
                kinetic_exp(space.dof{l},psi);
            end
            
            % Transition/permanent dipoles, polarizabilities
            if isfield(time,'pulse')
                
                % If possible, use field at next time step
                if k+1+time.steps.offset < length (time.efield.grid{1})
                    e = [ ...
                        time.efield.grid{1}(k+1+time.steps.offset) ...
                        time.efield.grid{2}(k+1+time.steps.offset)];
                else
                    e = [ ...
                        time.efield.grid{1}(k+time.steps.offset) ...
                        time.efield.grid{2}(k+time.steps.offset)];
                end
                pol_propa   (obj, psi, e);
                trans_propa (obj, psi, e, true);
                perm_propa  (obj, psi, e);
            end
            
            % Potential energy
            pot_propa (obj, psi);
            
        end
        
    end
    
end
