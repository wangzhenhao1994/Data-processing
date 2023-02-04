%--------------------------------------------------------------------------
%
% Sampling densities by bundles of trajectories, i.e. swarms of points 
% in phase space, the evolution of which is governed by classical dynamics. 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-.... Burkhard Schmidt
%
% see the README file for license details.

classdef traj < generic & handle
    
    properties (Access = public)
        
        n_p         % number of phase space points/particles
        
        pos         % positions of particles (cell vector)
        mom         % momenta of particles (cell vector)
        frc         % forces acting on "particles" (cell vector)
        
        pot         % potential energy
        kin         % kinetic energy
        
        pot_mat     % diabatic potential matrices
        frc_mat     % diabatic forces
        
        cha         % Assigning channel indices to each trajectory
        
        D_new       % Adiabatic energies (eigenvalues)
        U_new       % Adiabatic states (eigenvectors) 
        U_old       % Adiabatic states (eigenvectors) of the previous step
        
        rnd_seed    % Toggle use of a predictable sequence of random numbers
        
        psi         % Associating quantum state vectors with each trajectory
        psi_old     % "Old" quantum state vectors for each trajectory
        psi_new     % "New" quantum state vectors for each trajectory
        
        choice_eq_motion   % Which equations of motion should be propagated?
        extra_tdse_solving % Needs the el. TDSE to be solved additionally? (depends on the choice_eq_motion and method)
        mom_kin            % kinetic momenta of particles (cell vector)
        frc_qt_adi         % Forces governing dynamics for adiabatic quantum trajectories
        frc_qt_dia         % Forces governing dynamics for diabatic quantum trajectories
        
        % Six text strings: Details of trajectory methods
        string1     % Sampling by classical trajectories
        string2     % Number of trajectories
        string3     % Seeding of random number generator
        string4     % Surface hopping variants
        string5     % Gradient descent details (if available)
        string6     % Scaling velocities|momenta (if available)
        % Seventh string occurring in plots from propagators, see folder +tmp/+traj
        
    end
        
    methods (Access = public)
        
        % Constructor: Setting default values
        function obj = traj (n,seed)
          
            global info
            
            % Inherit from constructor of generic superclass
            obj = obj@generic;
            
            % Number of phase space points/particles
            if isempty(n)
                obj.n_p = 1000;          % number of phase space points/particles
            else
                obj.n_p = n;
            end
            
            % Reproducable sequence of random numbers
            if isempty(seed)
                obj.rnd_seed = [];        
            else
                obj.rnd_seed = seed;
                if strcmpi(info.system,'Matlab') 
                  rng(obj.rnd_seed);
                elseif strcmpi(info.system,'Octave')
                  randn("state",obj.rnd_seed) % initialize generation of normally distributed random numbers
                end
            end
            
            obj.choice_eq_motion = 'cl';    % usual classical dynamics
            obj.extra_tdse_solving = 0;     % Do not solve the el. TDSE if not needed
            
            % Text string => setting the name of this class
            obj.string0 = 'traj';
            
            % Text strings => lower right corner of expect-plots
            obj.string1 = 'Sampling densities by swarms of trajectories';
            obj.string2 = ['Number of particles : ' int2str(obj.n_p)];
            if isempty(obj.rnd_seed)
                obj.string3 = 'Arbitrary sequence of random numbers';
            else
                obj.string3 = ['Reproducable sequence of random numbers; seed : ' int2str(obj.rnd_seed)];
            end
            obj.string4 = 'Purely classical trajectories: No surface hopping';
            obj.string5 = ' ';
            obj.string6 = ' ';
        end
         
        % Logfile/console output
        function disp(obj)
            prt.disp (' ')
            prt.disp ('***************************************************************')
            prt.disp ('Solving the (quantum-)classical Liouville equation using ')
            prt.disp ('sampling of phase-space densities by trajectory ensembles')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp (obj.string1);
            prt.disp (obj.string2);
            prt.disp (obj.string3);
            prt.disp (' ')
            prt.disp (obj.string4);
            prt.disp (obj.string5);
            prt.disp (obj.string6);
            prt.disp (' ')
        end
        
        % How much memory needed to save a "traj" object
        function out = memory_size (obj)
            global space
            out = obj.n_p * ( space.n_dim * 16 + 4 );
        end        
        
        % Continuity of eigenvectors U: avoid sudden sign changes of columns
        function ensure_continuity_eig (obj, U_old)
            global hamilt
            
            n_traj = length(obj.U_new(1,1,:));
            
            for p = 1:hamilt.coupling.n_eqs
                minus = zeros(n_traj,1);
                minus(:) = sum(obj.U_new(:,p,:) .* U_old(:,p,:)) < 0;
                ind_minus = find(minus==1);
                
                if( isempty(ind_minus) == 0 )
                    obj.U_new(:,p,ind_minus) = - obj.U_new(:,p,ind_minus);
                end
            end
        end
        
        function save_previous(obj)
            obj.U_old = obj.U_new;
        end
        
        % see separate files for the following public methods
        init_obj  ( obj )                % Initial conditions
        propagate ( obj, step )          % Propagation
        observe   ( obj, step )          % Expectation values and uncertainties
        init_ham  ( obj )                % Initialization of Hamiltonian
        apply_ham ( obj )                % Application of Hamiltonian
        adiabatic ( obj, step, direction)% Adiabatic<=>diabatic transformation
        save      ( obj, step )          % Saving dynamics to file(s)
        load      ( obj, step )          % Loading dynamics from file(s)
        save_0    ( obj       )          % Saving general info to file
        load_0    ( obj, choice )        % Loading general from file
        
        % electronic TDSE propagators
        tdse_adi (obj, step_size, D, obj_U_new, obj_U_old)
        tdse_dia (obj)
        
    end
    
end

