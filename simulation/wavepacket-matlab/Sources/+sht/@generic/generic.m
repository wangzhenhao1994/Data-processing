%--------------------------------------------------------------------------
%
% Generic class for all "surface hopping trajectories" class definitions
%
% For use in mixed quantum-classical dynamics where the trajectories
% may undergo transitions between different (diabatic or adiabatic)
% states of the quantum system in a stochastic manner.
%
% Note that the method traj_hop defined here is calling methods prob_hop 
% and prep_hop which have to be implemented in each of the sub-classes.
% In addition, the sub-classes may overwrite the methods init_hop as well
% as prep_hop, thus providing the necessary flexibility in setting up 
% implementations of different variants of surface hopping.
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

classdef generic < traj & handle
    
    properties (Access = public)
        
        Allowed     % Indices of trajectories where hopping is allowed
        Forbidden   % Indices of trajectories where hopping is forbidden
        
        rescale     % Toggle momentum rescaling upon surface hopping
        sca_nac     % Scaling along non-adiabatic coupling vector
        
        verbose     % Toggle diagnostic extra output
        
        ind_h       % Indices of trajectories about to hop in current time step
        
        consec_hops % allow only consecutive hops
        
        ham         % Associating Hamiltonian matrices with each trajectory
                
        allowed     % Indices of trajectories where hopping is allowed
        forbidden   % Indices of trajectories where hopping is forbidden
        
    end
    
    methods (Access = public)
        
        % Constructor: Setting defaults and a text string
        function obj = generic (n,seed)
            obj = obj@traj(n,seed);     % Inherit from superclass
            obj.string0 = 'sht';        % Setting the name of this class
            obj.verbose = false;        % Diagnostic extra output
            obj.consec_hops = false;    % hops between all levels are allowed
            obj.string4 = 'Superclass for all surface hopping variants';
            obj.string5 = 'No gradient descent employed';
        end
                
        % Initialization of hopping
        function init_obj (obj)
            global hamilt
            
            init_obj@traj(obj);
            
            % Hamiltonian matrices attached to trajectories
            obj.ham       = cell(hamilt.coupling.n_eqs);
            for m=1:hamilt.coupling.n_eqs
                for n=1:hamilt.coupling.n_eqs
                    obj.ham {n,m} = zeros(obj.n_p,1);
                end
            end
            
            % Indices of trajectories where hopping is allowed/forbidden
            obj.Allowed   = cell(hamilt.coupling.n_eqs);
            obj.Forbidden = cell(hamilt.coupling.n_eqs);
                        
            % Toggle momentum rescaling upon surface hopping
            if isempty(obj.rescale)
                if strcmpi(hamilt.coupling.represent,'dia')
                    obj.rescale = false;
                elseif strcmpi(hamilt.coupling.represent,'adi')
                    obj.rescale  = true;
                end
            end
            
            % Scaling along non-adiabatic coupling vectors
            % Else, scaling is along the momentum vectors
            if obj.rescale
                if isempty(obj.sca_nac)
                    obj.sca_nac = false;
                end
            end
            
            % Output scaling method to text strings => plot of expectation values
            if obj.rescale
                if obj.sca_nac
                    obj.string6 = 'Rescaling momenta along NAC vectors';
                else
                    obj.string6 = 'Rescaling momenta along old momenta';
                end
            else
                obj.string6 = 'No rescaling of momenta';
            end
            
            % Scaling along NAC vectors only in adiabatic representation
            if obj.rescale
                if strcmpi(hamilt.coupling.represent,'dia') && obj.sca_nac
                    prt.error ('Scaling along NAC vectors only in adiabatic representation')
                end
            end
            
        end
        
        % Main class for performing surface hopping
        function hop_main( obj, first_call )
            
            prep_hop (obj, first_call)       % Preprocessing
            traj_hop (obj)                   % Trajectory hopping
            
        end
        
        % Preprocessing: before hopping
        function prep_hop ( obj,first_call )
            global hamilt
            
            % Reset statistics for allowed/forbidden hopping events
            % So far, this is used in scatter plots only
            if first_call
                obj.Allowed   = cell (hamilt.coupling.n_eqs);
                obj.Forbidden = cell (hamilt.coupling.n_eqs);
            end
            
            % Update Hamiltonians attached to Q/C trajectories
            eval_ham ( obj )
            
        end

         
        function store_allowed_forbidden (obj, allowed, forbidden, m, n)
            % Unique values in array
            obj.Allowed  {n,m} = union ( obj.Allowed  {n,m} , allowed   );
            obj.Forbidden{n,m} = union ( obj.Forbidden{n,m} , forbidden );
        end
        
        
        % see separate files for the following public methods
        eval_ham  ( obj, nac )           % Evaluate Hamiltonian operators
        propagate ( obj, step )          % Propagation
        traj_hop  ( obj )                % Perform the trajectory hopping
        after_hop ( obj, mom_new,ind_h_allowed,allowed,~,n)
        disp_hop  ( obj, pos_cell,mom_cell,mom_new_cell,D,ind_h_allowed,ind_h_forbidden,allowed,m,n)
        [pos,mom,pot_mats,frc_mats,U,D] = get_quantities    ( obj, ind,~,~)
        [mom_new, allowed, forbidden]   = mom_rescaling     ( obj, mom,pot_mats,frc_mats,U,D,m,n)
        [mom_new, allowed, forbidden]   = mom_rescaling_vec ( obj, mom,D,vec,m,n)
        
    end
        
end
    
