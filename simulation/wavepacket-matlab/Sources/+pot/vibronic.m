%--------------------------------------------------------------------------
% Potential energy matrix for vibronic coupling 
%--------------------------------------------------------------------------
% 
% For a one-dimensional crystal (i.e. a chain with periodic boundaries) 
% or, equivalently, for a regular polygon-shaped molecule (C_n symmetry)        
% Description of excitons in terms of a Frenkel site Hamiltonian
% with neutral sites assumed to have two orbitals (HOMO and LUMO)
% Including electron-phonon coupling (Jahn-Teller and pseudo-Jahn-Teller) 
% within linear approximation.
%
% A. Viel and W. Eisfeld, J. Chem. Phys. 122, 204317 (2005)     
% doi:10.1063/1.1904594
%
% F. Giustino, Rev. Mod. Phys. 89 (1), 1-63 (2017)
% doi:10.1103/RevModPhys.89.015003

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2018 Burkhard Schmidt's group
%               2009 Burkhard Schmidt
%
% see the README file for license details.

classdef vibronic < pot.generic & handle
    
    properties (Access = public)
        
        N           % Number of electronic sites
        U           % Electronic site energy 
        W           % Electronic coupling strength
                
        O           % Harmonic frequency (omega) of one particle
        S           % Selection of phonons
        
        K           % Linear electron-phonon coupling strength ( JT)
        L           % Linear electron-phonon coupling strength (PJT)
       
    end
        
    properties (Access = private)
        
        e           % Electronic energies
           
        j           % Wavenumber indices
        k           % Wavenumbers in units of 1/a
        l           % Labels with wavenumbers

        o           % Phonon harmonic frequencies
        r           % Frequency ratio
        
        q           % Linear electron-phonon coupling
        
        v1          % First Taylor coefficients
        v2          % Second Taylor coefficients (force constants)

    end
    
    methods (Access = public)
        
       % Constructor: Set default values
        function obj = vibronic
            
            % Default: consider all(!) phonon E modes
            obj.N = 3;                   % prototypical case: (E+A) x (e+a)
            obj.U = 0;                   % without restriction of generality
            obj.K = 0;                   % no JT coupling
            obj.L = 0;                   % no PJT coupling
            obj.S = 1:obj.N;             % consider ALL phonon modes
            
        end
        
        
        % Initialize potential: Set/check parameters
        function init (obj)
            global hamilt space
            
            % Check number of electronic sites
            if hamilt.coupling.n_eqs ~= obj.N
                prt.error ('Number of electronic sites inconsistent')
            end
            
            % Number of doublets (E irrep)
            if mod(obj.N,2)              % If N is odd
                n = +(obj.N-1)/2;                      
            else                         % If N is even
                n = +(obj.N-2)/2;   
            end
            
            % Bloch states with discrete quasi-momentum/wavenumber k = 2*pi*j / (a*N)
            obj.j = zeros (obj.N,1);                          % Integer values
            obj.l = cell  (obj.N,1);                          % Labels (strings)
            obj.j (1) = 0;                                    % A irrep: k = 0
            obj.l {1} = 'a';                                  % A irrep
            for d=1:n
                obj.j (2*d  ) = +d;                           % E irreps
                obj.l {2*d  } = ['x' int2str(d)];       
                obj.j (2*d+1) = -d;                           % E irreps
                obj.l {2*d+1} = ['y' int2str(d)];       
            end
            if mod(obj.N,2) == 0
                obj.j(end)= n + 1;                            % B irrep: k = pi/a
                obj.l{end}= 'b';
            end
            obj.k = 2*pi*obj.j/obj.N;                         % Wavenumber (for lattice constant a=1)
            
            % Electronic energies and labels
            obj.e = obj.U + 2*obj.W*cos(obj.k);               % Eigenenergy of electron
            for m = 1:obj.N
                hamilt.coupling.labels{m} = lower(obj.l{m});
            end
            
            % Check selection of phonon coordinates
            if length (obj.S) ~= space.n_dim
                prt.error ('Number of selected phonons inconsistent')
            end
            if any ( obj.S < 1 )
                prt.error ('Phonon selection index too small')
            end
            if any ( obj.S > obj.N )
                prt.error ('Phonon selection index too large')
            end
            
            % Labeling the phonon coordinates
            for d = 1:length(obj.S)
                space.dof{d}.label = upper(obj.l{obj.S(d)});
            end
                
            % Phononic Bloch states with quasi-momentum/wavenumber q_i = 2*pi*ii / (a*N)
            obj.r = sqrt(2-2*cos(obj.k));                     % Ratio of frequencies of phonons
            obj.o = obj.O .* obj.r;                           % Harmonic frequencies of phonons
            
            % Linear electron-phonon couplings in complex (+/-) representation
%             C = cell (obj.N,1);
%             C{1} = eye(obj.N);                                % A irrep
%             for e=1:n
%                 C{2*e  } = circshift(C{1},+e,2);              % E irreps; Q_plus               
%                 C{2*e+1} = circshift(C{1},-e,2);              % E_irreps: Q_minus
%             end
%             if mod(obj.N,2) == 0
%                 C{end}= circshift(C{1},n+1,2);                % B irrep
%             end

            c = cell (obj.N,1);
            for d = 1:obj.N
                c{d} = zeros(obj.N);
            end
            for row = 1:obj.N
                jr = obj.j(row);
                for col = 1:obj.N
                    jc = obj.j(col);
                    if mod(jr-jc,obj.N)==0                              
                        c{1}(row,col) = 1;                    % A irrep
                    end
                    for d = 1:n                               % E irreps
                        if mod(jr-jc+d,obj.N)==0                          
                            c{2*d  }(row,col) = 1;            % Q_plus
                        end
                        if mod(jr-jc-d,obj.N)==0                     
                            c{2*d+1}(row,col) = 1;            % Q_minus
                        end
                    end
                    if mod(obj.N,2) == 0                      % B irrep
                        if mod(jr-jc+(n+1),obj.N)==0  
                            c{end}(row,col) = 1;        
                        end
                    end
                end
            end
            
            % Transform electrons from complex (+/-) to real (x/y) representation
            s = [+1 +1; +1i -1i]/sqrt(2);
            t = zeros(obj.N);
            t (1,1) = 1;                                      % A irrep
            for d = 1:n
                t(2*d:2*d+1, 2*d:2*d+1) = s;                  % E irreps
            end
            if mod(obj.N,2) == 0
                t(end,end) = 1;                               % B irrep
            end
            
            for d = 1:obj.N
                c{d} = t * c{d} * t';
            end
            
            % Transform phonons from complex (+/-) to real (x/y) representation
            obj.q = cell (obj.N,1);
            obj.q{1} = c{1};                                  % A irrep
            for d = 1:n
                obj.q{2*d  } = (c{2*d}+c{2*d+1}) ;            % E irreps: real part
                obj.q{2*d+1} = (c{2*d}-c{2*d+1})/1i;          % E irreps: imag part 
            end
            if mod(obj.N,2) == 0
                obj.q{end} = c{end};                          % B irrep
            end

            % Getting Taylor coefficients for use in methods V and F (see below)
            obj.v2 = obj.o(obj.S).^2;                         % Force constants (assuming m=1)
            obj.v1 = cell (obj.N);
            for m = 1:obj.N
                for n = m:obj.N
                    obj.v1{m,n} = zeros(1,length(obj.S));     % Row vector: linear couplings
                    if abs(obj.j(m))==abs(obj.j(n))
                        f = obj.K;                            % Jahn-Teller
                    else
                        f = obj.L;                            % pseudo-Jahn-Teller
                    end
                    for d = 1:length(obj.S)
                        obj.v1{m,n}(d) = f * sign(math.real(obj.q{obj.S(d)}(m,n)));
                    end
                end
            end
            
        end

        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            if obj.row==1 && obj.col==1
                prt.disp ('Coupled electron-phonon dynamics: Jahn-Teller and beyond')
                prt.disp ('***************************************************************')
                prt.disp (' ')
                prt.disp ('**************')
                prt.disp ('Electrons only')
                prt.disp ('**************')
                prt.disp (' ')
                prt.disp (['Number of electronic sites    : ' int2str(obj.N)])
                prt.disp (['Electronic site energy        : ' num2str(obj.U)])
                prt.disp (['Electronic coupling strength  : ' num2str(obj.W)])
                prt.disp (' ')
                for d=1:obj.N
                    prt.disp (['Exciton wavenumber index      : ' int2str(obj.j(d))])
                    prt.disp (['Exciton wavenumber (for a=1)  : ' num2str(obj.k(d))])
                    prt.disp (['Exciton energy                : ' num2str(obj.e(d))])
                    prt.disp (' ')
                end
                prt.disp ('************')
                prt.disp ('Phonons only')
                prt.disp ('************')
                prt.disp (' ')
                prt.disp (['1-particle harmonic frequency : ' num2str(obj.O)])
                prt.disp (' ')
                for d=1:obj.N
                    prt.disp (['Phonon wavenumber index       : ' int2str(obj.j(d))])
                    prt.disp (['Phonon wavenumber (for a=1)   : ' num2str(obj.k(d))])
                    prt.disp (['Phonon harmonic frequency     : ' num2str(obj.o(d))])
                    prt.disp (['Frequency ratio               : ' num2str(obj.r(d))])
                    prt.disp (' ')
                end
                prt.disp ('Selection of phonons')
                prt.disp (' ')
                for d=1:length(obj.S)
                    prt.disp ([int2str(obj.S(d)) ': ' obj.l{obj.S(d)}])
                end
                prt.disp (' ')
                prt.disp ('************************')
                prt.disp ('Electron-phonon coupling')
                prt.disp ('************************')
                prt.disp (' ')
                prt.disp (['Coupling constant ( JT) : ' num2str(obj.K)])
                prt.disp (['Coupling constant (PJT) : ' num2str(obj.L)])
                prt.disp (' ')
                for d = 1:length(obj.S)
                    disp(strcat('Q_',obj.l(obj.S(d))))
                    disp(       real(obj.q{obj.S(d)}))
                end
            else
                prt.disp ('Same as above')
                prt.disp ('***************************************************************')
            end
            prt.disp (' ')
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)
            
            % Electronic energies
            if obj.row==obj.col
                V = obj.e(obj.row) * ones(size(r{1}));      
            else
                V = zeros(size(r{1}));
            end
            
            % Phonon vibrations
            if obj.row==obj.col
                for d=1:length(obj.S)
                    V = V + obj.v2(d) * r{d}.^2 / 2;
                end
            end
            
            % Linear electron-phonon-coupling
            for d=1:length(obj.S)
                V = V + obj.v1{obj.row,obj.col}(d) * r{d};
            end
            
        end

        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            
            F = cell(size(r));
            
            % Initialize forces
            for d=1:length(obj.S)
                F{d} = zeros(size(r{1}));
            end
            
            % Phonon vibrations
            if obj.row==obj.col
                for d=1:length(obj.S)
                    F{d} = F{d} - obj.v2(d) * r{d};
                end
            end
            
            % Linear electron-phonon-coupling
            for d=1:length(obj.S)
                F{d} = F{d} - obj.v1{obj.row,obj.col}(d);
            end

            
        end
        
    end

end