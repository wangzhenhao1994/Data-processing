%--------------------------------------------------------------------------
%
% Set up matrix representation of the Hamiltonian operator 
% Solve the corresponding eigenproblem, i.e. the 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef eigen < handle
    
    properties (Access=public)
        
        start       % Index of first eigenstate wanted
        stop        % Index of last eigenstate wanted
        number      % Number of eigenstates wanted
        storage     % Full or sparse storage of Hamiltonian matrix
        cutoff      % Set matrix entries below this threshold to zero
        symmetry    % Symmetry (irr. rep.) label
        matrix      % Matrix representation of Hamiltonian
        transform   % Transformation matrix (if symmetry adation)
        
        eig_vals    % Eigenvalues
        eig_vecs    % Eigenvector
        
    end
    
    methods (Access=public)
        
        % Constructor: Set default values
        function obj = eigen

            obj.start   = 00;  % ground state
            obj.stop    = 10;  % 10th excited state
            obj.cutoff  = 0.0; % no cutoff
            obj.storage = 'f'; % full
            
        end
        
        function init (obj)
            global time
            
            % Artificial mapping of eigenstates onto time step 
            time.steps.m_start  = obj.start;
            time.steps.m_stop   = obj.stop;
            obj.number          = obj.stop - obj.start + 1;
            time.steps.m_number = obj.number;
            time.steps.m_grid   = (obj.start : obj.stop)';
            time.steps.t_total  = obj.stop;
            time.steps.s_delta  = 1e-10; % dummy setup to avoid crashes

        end
            
        
        % Display objects of this class
        function disp (obj)
            global hamilt space
            
            prt.disp('***************************************************************')
            prt.disp('Solve TISE by diagonalization of the Hamiltonian  ')
            prt.disp('***************************************************************')
            prt.disp('  ')
            prt.disp(['Size of matrix                : ' int2str(space.n_tot * hamilt.coupling.n_eqs) '^2'])
            if obj.cutoff == 0
                prt.disp('No cutoff value');
            else
                prt.disp(['Cut-off absolute values below : ' num2str(obj.cutoff)])
            end
            if obj.storage == 'f'
                prt.disp('Storing full matrix');
            elseif obj.storage == 's'
                prt.disp('Storing sparse matrix');
            else
                prt.error('Wrong storage scheme. hamilt.eigen.storage has to be "f" (full) or "s" (sparse)');
            end
            prt.disp('   ')
            prt.disp(['First eigenstate investigated : ' num2str(obj.start)])
            prt.disp(['Last  eigenstate investigated : ' num2str(obj.stop)])
            
            
        end
        
        setup (obj)
        symm (obj)
        diag (obj)
        
    end
    
end

