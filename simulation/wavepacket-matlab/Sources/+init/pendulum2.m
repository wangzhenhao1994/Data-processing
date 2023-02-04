%------------------------------------------------------------------------------
%
% Stationary states of a generalized planar pendulum
%
%              V(R) = - eta*cos(R) - zeta*cos(R)^2
% 
% This class gives the eigenfunctions of a generalized planar pendulum
% the time-independent Schrödinger equation of which classifies as
% quasi-exactly solvable and conditionally exactly solvable, i.e. only
% a certain number of eigenstates can be found for integer values of kappa
% yielding parameters eta = kappa*beta and zeta = beta^2. Note that the 
% potential parameters can be either given as eta,zeta or kappa,beta.
% For more complete explanations, see our work in the following references: 
%
% S. Becker, M. Mirahmadi, B. Schmidt, K. Schatz, B. Friedrich 
% Eur. Phys. J. D, 71(6), 149 (2017)
% https://dx.doi.org/10.1140/epjd/e2017-80134-6
%
% B. Schmidt, & B. Friedrich
% Front. Phys., 2, 1–16 (2014)
% https://dx.doi.org/10.3389/fphy.2014.00037
%
% The irreducible representations A1, A2 give 2pi-periodic solutions with
% even and odd parity, respectively. The irreducible representations B1, B2 
% give 2pi-antiperiodic solutions with even and odd parity, where the value 
% of parity will change depends on the sign of beta. Note that to get the 
% exact eigenvalues for irreps. B1 and B2 you should use a 4*pi interval 
% for angle R whereas a 2*pi interval suffices for A1 and A2 states.
%
% Because this problem is only conditionally quasi-exactly solvable, 
% the number of analytical solutions is limited. Our built-in error 
% messages show the restrictions on available analytical states with 
% quantum number n_q.
%
% This classdef can also create a von-Mises distribution exp(beta cos X) 
% for special values: n_q = 0, kappa = 1, irrep = 'A1' .
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019 Burkhard Schmidt & Marjan Mirahmadi
%
% see the README file for license details.

classdef pendulum2 < init.generic & handle
    
    properties (Access = public)
        
        beta    % Width parameter
        kappa   % number of solutions
        eta     % Potential parameter
        zeta    % Potential parameter
        n_q     % Index of analytic eigenstate (starts from zero)
        irrep   % irreducible representation A1 A2 B1 B2
        
    end
    
    properties (Access = private)
        
        V         % Eigenvectors of matrix representation of pendular Hamiltonian
        D         % Eigenvalues of matrix representation of pendular Hamiltonian
        
    end
    
    methods (Access = public)
        % Constructor: Set default values
        function obj = pendulum2
            
        end
        
        % Initialize GP functions: Set/check parameters
        function init (obj)
            global space
            
            % FFT only: get the coordinate interval 
            if ~isa(space.dof{obj.dof},'dof.fft')
                prt.error ('This wavefunction allowed only for FFT grids!')
            end
            dx = space.dof{obj.dof}.x_max - space.dof{obj.dof}.x_min;
            
            % Get beta and kappa if not provided by the user
            if isempty(obj.beta) || isempty(obj.kappa)
                obj.beta =  sqrt(obj.zeta);
                obj.kappa = obj.eta / obj.beta;
            end
            
            % Check for integer kappa
            if round(obj.kappa)~=obj.kappa
                    prt.error ('kappa should be integer')
            end
            
            if (obj.kappa<1)
                prt.error ('kappa should be positive')
            end
            
                       
            switch obj.irrep
                case 'A1'
                    if mod(dx,2*pi)
                        prt.error ('A1 irrep available only for 2pi interval')
                    end
                    if ~mod(obj.kappa,2)
                        prt.error ('A1 irrep is only for odd kappa')
                    end
                    T = A1(obj); % Setup matrix T
                    
                case 'A2'
                    if mod(dx,2*pi)
                        prt.error ('A2 irrep available only for 2pi interval')
                    end
                    if ~mod(obj.kappa,2)
                        prt.error ('A2 irrep is only for odd kappa')
                    end
                    if obj.kappa < 3
                        prt.error ('A2 irrep not for kappa=1')
                    end
                    T = A2(obj); % Setup matrix T
                    
                case 'B1'
                    if mod(dx,4*pi)
                        prt.error ('B1 irrep available only for 4pi interval')
                    end
                    if mod(obj.kappa,2)
                        prt.error ('B1 irrep is only for even kappa')
                    end
                    T = B1(obj); % Setup matrix T
                    
                case 'B2'
                    if mod(dx,4*pi)
                        prt.error ('B2 irrep available only for 4pi interval')
                    end
                    if mod(obj.kappa,2)
                        prt.error ('B2 irrep is only for even kappa')
                    end
                    T = B2(obj); % Setup matrix T
                    
            end
            
            % Eigenvectors (V) and eigenvalues (D) of matrix T
            [obj.V,obj.D] = eig (T,'vector');
            if obj.beta>0
                obj.D = flipud(obj.D);
                obj.V = fliplr(obj.V);
            end
            
            % Check whether "index" is inside analytic part of spectrum
            if obj.n_q > length (T)-1
                prt.error ('The quantum number is too large for this kappa!')
            end
            
        end
        
        % Display GP function, overloading default disp method
        function disp(obj)
            
            prt.disp ('Analytic eigenfunctions of a generalized planar pendulum                                ')
            prt.disp ('***************************************************************')
            prt.disp (' ' )
            prt.disp ( ['Potential parameter beta   : ' num2str(obj.beta)] )
            prt.disp ( ['Potential parameter kappa  : ' int2str(obj.kappa)] )
            prt.disp (' ' )
            prt.disp ( ['Potential parameter eta    : ' num2str(obj.beta*obj.kappa)] )
            prt.disp ( ['Potential parameter zeta   : ' num2str(obj.beta^2)] )
            prt.disp (' ' )
            prt.disp ( ['Irreducible representation : '         obj.irrep] )
            prt.disp (' ' )
            prt.disp ( ['Quantum number             : ' int2str(obj.n_q)] )
            prt.disp ( ['Eigenenergy                : ' num2str(-obj.D(obj.n_q+1))])
            
        end
        
        % Evaluate wave function on a grid  ==> Q/M propagation
        function wave (obj)
            global space
            
            % Adding monomials
            obj.dvr = zeros (size(space.dvr{obj.dof}));
            for j = 0:length (obj.D)-1
                obj.dvr = obj.dvr + obj.V(j+1,obj.n_q+1) .* cos(1/2*space.dvr{obj.dof}).^(2*j);
            end
            
            % Multiply with seed function imposing the right symmetry
            switch obj.irrep
                case 'A2'
                    obj.dvr = obj.dvr .* sin(space.dvr{obj.dof});
                case 'B1'
                    obj.dvr = obj.dvr .* cos(space.dvr{obj.dof}/2);
                case 'B2'
                    obj.dvr = obj.dvr .* sin(space.dvr{obj.dof}/2);
            end
            
            % Multiply with von Mises function
            obj.dvr = obj.dvr .* exp(obj.beta * cos(space.dvr{obj.dof}));
        end
        
        % Sample phase space density ==> Q/C propagation
        function traj (~, ~)
            prt.error ('Code for phase space sampling still missing')
        end
        
        
    end
    
    methods (Access = private)
        
        function T = A1 ( obj )
            k = obj.kappa;
            b = obj.beta;
            N = (k-1)/2;
            ell = 0:N;   T =     diag (b^2 - ell.^2 + 4*b*ell - (k-1)*b);     % diagonal
            ell = 1:N;   T = T + diag (ell.^2 - ell/2, +1);                   % super-diagonal
            ell = 0:N-1; T = T + diag (-4*b*ell + 2*(k-1)*b, -1);              % sub-diagonal
        end
        
        function T = A2 ( obj )
            k = obj.kappa;
            b = obj.beta;
            N = (k-3)/2;
            ell = 0:N;   T =     diag (b^2 - ell.^2 + 4*b*ell -(k-3)*b - 2*ell - 1); % diagonal
            ell = 1:N;   T = T + diag (ell.^2 + ell/2, +1);                 % super-diagonal
            ell = 0:N-1; T = T + diag (-4*b*ell + 2*(k-3)*b, -1);            % sub-diagonal
        end
        
        function T = B1 ( obj )
            k = obj.kappa;
            b = obj.beta;
            N=k/2-1;
            ell = 0:N;   T =     diag (b^2 - ell.^2 + 4*b*ell - (k-3)*b - ell - 1/4); % diagonal
            ell = 1:N;   T = T + diag (ell.^2 + ell/2, +1);                % super-diagonal
            ell = 0:N-1; T = T + diag (-4*b*ell + 2*(k-2)*b, -1);           % sub-diagonal
        end
        
        function T = B2 ( obj )
            k = obj.kappa;
            b = obj.beta;
            N=k/2-1;
            ell = 0:N;   T =     diag (b^2 - ell.^2 + 4*b*ell - (k-1)*b - ell - 1/4); % diagonal
            ell = 1:N;   T = T + diag (ell.^2 - ell/2, +1);                  % super-diagonal
            ell = 0:N-1; T = T + diag (-4*b*ell + 2*(k-2)*b, -1);             % sub-diagonal
        end
        
    end
    
end