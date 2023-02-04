%--------------------------------------------------------------------------
% Homogenous chain of harmonic oscillators (model for phonons in 1D)
%
% For analytic solutions of the time independent Schroedinger equation, see
%     https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator
%     https://doi.org/10.1063/5.0074948 (Gelss,  Klein, Matera, Schmidt)
%
% For analytic solutions of the time dependent Schroedinger equation, see
%     https://doi.org/10.1016/0370-1573(92)90093-F (Scott: Davydov soliton)
%
% Can be used also for coupled electrons and phonons where
% tuning (chi,rho,sig) and coupling (tau) mechanisms involve
% nearest-neighbor interactions only
%
% optionally with periodic boundary conditions (chain==>ring)
% optionally with position restraints
% optionally with center-of-mass restraints
% optionally with electron-phonon-tuning|coupling
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019 Burkhard Schmidt
%
% see the README file for license details.

classdef chain < pot.generic & handle
    
    properties (Access = public)
        
        oNN         % harmonic frequency: nearest neighbors
        oPR         % harmonic frequency: position restraints
        oCM         % harmonic frequency: center-of-mass restraint
        r_e         % equilibrium distance
        pbc         % toggle periodic boundary conditions
        off         % vertical energy offset
        chi         % electron-phonon-tuning: localized
        rho         % electron-phonon-tuning: non-symmetric
        sig         % electron-phonon-tuning: symmetrized
        tau         % electron-phonon-coupling: pair distance
                
    end
    
    methods (Access = public)
        
       % Constructor: Set default values
        function obj = chain
            obj.oNN = 0;
            obj.oPR = 0;
            obj.oCM = 0;
            obj.r_e = 0;
            obj.pbc = false;
            obj.off = 0;
            obj.chi = 0;
            obj.rho = 0;
            obj.sig = 0;
            obj.tau = 0;
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            global space
            if obj.oNN < 0
                prt.error ('Invalid harmonic frequency: oNN')
            end
            if obj.oPR < 0
                prt.error ('Invalid harmonic frequency: oPR')
            end
            if obj.oCM < 0
                prt.error ('Invalid harmonic frequency: oCM')
            end
            if obj.r_e < 0
                prt.error ('Invalid equilibrium distance')
            end
            if ~islogical (obj.pbc)
                prt.error ('Invalid boolean: periodic boundary conditions')
            end
            if space.n_dim<2
                prt.error ('At least two DoFs required')
            end
            for k=2:space.n_dim
                if space.dof{k}.n_pts ~= space.dof{k-1}.n_pts
                    prt.error ('Equal number of grid points for every DoF required')
                end
            end
            for k=2:space.n_dim
                if space.dof{k}.mass ~= space.dof{k-1}.mass
                    prt.error ('Equal mass for every particle required')
                end
            end
            if obj.sig && space.n_dim<3
                prt.error('Electron-phonon-tuning variant 2 only for more than 2 sites')
            end
            
            if obj.chi~=0 || obj.rho~=0 || obj.sig~=0 
               if obj.row ~= obj.col
                   prt.error('electron-phonon tuning for diagonal only')
               end
            end

            if obj.tau~=0 
               if ~( obj.row+1==obj.col || (obj.row==1&&obj.col==space.n_dim) )
                   prt.error('electron-phonon coupling for super-diagonal or upper right only')
               end
            end

        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            global space
            disp@pot.generic (obj)
            prt.disp ('Chain of harmonic oscillators, optionally with PBC')
            prt.disp ('Optionally with electron-phonon-tuning|coupling')
            prt.disp ('***************************************************************')
            prt.disp ( ' ' )
            prt.disp ( [ 'Number of sites in chain|ring             : ' int2str(space.n_dim) ])
            prt.disp ( [ 'Mass of particles (from grid parameters)  : ' num2str(space.dof{1}.mass) ] )
            prt.disp ( [ 'Harmonic vibration frequency NN           : ' num2str(obj.oNN) ] )
            prt.disp ( [ 'Harmonic vibration frequency PR           : ' num2str(obj.oPR) ] )
            prt.disp ( [ 'Harmonic vibration frequency CM           : ' num2str(obj.oCM) ] )
            prt.disp ( [ 'Equilibrium chain distance                : ' num2str(obj.r_e) ] )
            prt.disp ( [ 'Periodic boundary conditions              : ' int2str(obj.pbc) ] )
            prt.disp ( [ 'Vertical energy offset (alpha/beta/eta)   : ' num2str(obj.off) ] )
            prt.disp ( [ 'electron-phonon-tuning: localized   (chi) : ' num2str(obj.chi) ] )
            prt.disp ( [ 'electron-phonon-tuning: non-symmetr (rho) : ' num2str(obj.rho) ] )
            prt.disp ( [ 'electron-phonon-tuning: symmetrized (sig) : ' num2str(obj.sig) ] )
            prt.disp ( [ 'electron-phonon-coupling: pair dist.(tau) : ' num2str(obj.tau) ] )
            prt.disp ( ' ' )
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)
            global space
            
            % Number of sites / DoF's
            N = length(r);
            
            % Summing up contributions from each component of position vector
            V = zeros(size(r{1}));
            
            % Nearest neighbor interactions
            if obj.oNN
                frc_const = space.dof{1}.mass/2 * obj.oNN^2; % force constant
                for k = 1:N-1
                    V = V + frc_const * ( r{k+1}-r{k} - obj.r_e ).^2 / 2;
                end
                
                if obj.pbc % Periodic boundary conditions
                    V = V + frc_const * ( r{N}-r{1} - (N-1) * obj.r_e ).^2 / 2;
                end
            end
            
            % Position restraints
            if obj.oPR 
                frc_const = space.dof{1}.mass * obj.oPR^2; % force constant
                for k = 1:N
                    V = V + frc_const * ( r{k} - k * obj.r_e ).^2 / 2;
                end
            end
            
            % Restraining the center of mass
            if obj.oCM 
                frc_const = N * space.dof{1}.mass * obj.oCM^2; % force constant
                CoM = zeros(n);
                for k = 1:N
                    CoM = CoM + r{k};
                end
                CoM = CoM / N;
                V = V + frc_const * ( CoM ).^2 / 2;
            end
            
            % Vertical energy offset (alpha or beta)
            if obj.off
                V = V + obj.off;
            end
            
            % Electron-phonon-tuning variant 0
            % Linear in local site coordinate
            if obj.chi
                k = obj.row;
                V = V + obj.chi * ( r{k} - k * obj.r_e );      
            end
            
            % Electron-phonon-tuning variant 1
            % Linear in nearest-neighbor distance
            if obj.rho
                k = obj.row;
                if k<N
                    V = V + obj.rho * ( r{k+1}-r{k} - 1 * obj.r_e );
                else
                    if obj.pbc % Periodic boundary conditions
                        V = V + obj.rho * ( r{1}-r{N} + (N-1) * obj.r_e );
                    end
                end
            end

            % Electron-phonon-tuning variant 2
            % Linear in second-nearest-neighbor distance
            if obj.sig 
                k = obj.row; 
                if k>1 && k<N
                    V = V + obj.sig * ( r{k+1}-r{k-1} - 2 * obj.r_e );
                elseif k==1
                    if obj.pbc % Periodic boundary conditions
                        V = V + obj.sig * ( r{2}-r{N} + (N-2) * obj.r_e );
                    end
                elseif k==N
                    if obj.pbc % Periodic boundary conditions
                        V = V + obj.sig * ( r{1}-r{N-1} + (N-2) * obj.r_e );
                    end
                end
            end
        
            % Electron-phonon-coupling variant 1
            % Linear in nearest-neighbor distance
            if obj.tau
                k1 = obj.row;
                k2 = obj.col;
                % For N=2 with PBC, the two if-blocks will cancel each other
                if k2-k1==1 % nearest neighbor
                    V = V + obj.tau * ( r{k2}-r{k1} -   1   * obj.r_e );
                end
                if obj.pbc && k1==1 && k2==N % periodic boundaries
                    V = V + obj.tau * ( r{1}-r{N}   + (N-1) * obj.r_e );
                end
            end

        end
            
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            global space
            
            if obj.r_e ~=0
                prt.error('Code missing here')
            end
            
            % Allocate memory
            N = length(r);
            F = cell(size(r));
            
            % Initialize forces
            for d=1:N
                F{d} = zeros(size(r{1}));
            end
            
            % Nearest neighbor interactions
            if obj.oNN
                frc_const = space.dof{1}.mass/2 * obj.oNN^2; % force constant
                
                % First site
                F{1} = F{1} + frc_const * ( r{2} - r{1} );
                if obj.pbc % Periodic boundary conditions
                    F{1} = F{1} - frc_const * ( r{1}-r{N} );
                end
                
                % Inner sites
                for k = 2:N-1
                    F{k} = F{k} + frc_const * ( r{k+1} - 2*r{k} + r{k-1} );
                end
                
                % Last site
                F{N} = F{N} - frc_const * ( r{N} - r{N-1});  
                if obj.pbc % Periodic boundary conditions
                    F{N} = F{N} + frc_const * ( r{1}-r{N} );
                end
            end
            
            % Position restraints 
            if obj.oPR
                frc_const = space.dof{1}.mass * obj.oPR^2; % force constant
                for k = 1:N
                    F{k} = F{k} - frc_const * r{k};
                end
            end
            
            % Restraining the center of mass
            if obj.oCM
                frc_const = N * space.dof{1}.mass * obj.oCM^2; % force constant
                for k = 1:N
                    F{k} = F{k} - frc_const/N * r{k};
                end
            end
            
            
            % Electron-phonon-tuning variant 0
            % Linear in local site coordinate
            if obj.chi
                k = obj.row;
                F{k} = F{k} - obj.chi;
            end
            
            % Electron-phonon-tuning variant 1
            % Linear in nearest-neighbor distance
            if obj.rho
                k = obj.row;
                if k<N
                    F{k}   = F{k}   + obj.rho;
                    F{k+1} = F{k+1} - obj.rho;
                else
                    if obj.pbc % Periodic boundary conditions
                        F{N}   = F{N} + obj.rho;
                        F{1}   = F{1} - obj.rho;
                    end
                end
            end
            
            % Electron-phonon-tuning variant 2
            % Linear in second-nearest-neighbor distance
            if obj.sig
                if k>1 && k<N
                    k = obj.row;
                    F{k-1} = F{k-1} + obj.sig;
                    F{k+1} = F{k+1} - obj.sig;
                elseif k==1
                    if obj.pbc % Periodic boundary conditions
                        F{N} = F{N} + obj.sig;
                        F{2} = F{2} - obj.sig;
                    end
                elseif k==N
                    if obj.pbc % Periodic boundary conditions
                        F{N-1}   = F{N-1} + obj.sig;
                        F{  1}   = F{  1} - obj.sig;
                    end
                end
            end
            
            % Electron-phonon-coupling variant 1
            % Linear in nearest-neighbor distance
            if obj.tau
                k1 = obj.row;
                k2 = obj.col;
                if k2-k1==1 % nearest neighbor 
                    F{k1} = F{k1} + obj.tau;
                    F{k2} = F{k2} - obj.tau;
                end
                if obj.pbc && k1==1 && k2==N % periodic boundaries
                    F{N} = F{N} + obj.tau;
                    F{1} = F{1} - obj.tau;
                end
            end
        end
    end
end