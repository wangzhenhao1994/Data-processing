%----------------------------------------------------------------------
%
% Real time Chebychev propagator for a quantum state vector
% ===========================================================
%
% Epansion of the real time evolution operator exp (-i H t / hbar )
% in imaginary Chebychev polynomials for time-independent Hamiltonian.
%
% R. Kosloff, J. Phys. Chem. 92(8), 2087, (1988)
%
%                    (    de         )        (              )
% exp(-i*H*dt) = exp (-i*(--+e   )*dt)  * exp (-i*alpha*H    )
%                    (    2   min    )        (          norm)
%
% where
%         2    (       de       )
% H     = -- * (H - I*(--+e   ) )
%  norm = de   (       2   min  ) 
%
%          de * dt
% alpha  = -------        
%             2
%
% DT      time step
% EMIN    minimum of the Hamiltonian
% DE      range of the Hamiltonian
% ALPHA   dimensionless parameter 
% Hnorm   Hamiltonian normalized to [-1,1]
% H       original Hamiltonian 
% I       unity operator
%
% Then the exponential of the time evolution operator 
% is expanded in a series of complex Chebychev polynomials.
%
%     (               )     N
% exp ( -i*alpha*H    )  = Sum a (alpha) * phi (-i*H    )
%     (           norm)    n=0  n             n     norm 
%  
% where the coefficients a_n are Bessel functions and
% where the phi_n are the complex Chebychev polynomials
% which are calculated using the recursion 
%
%    phi  = I   (unity)
%       0
%
%    phi  = -i*H
%       1       norm
%
%    phi  = -2*i*H     phi    + phi    , n>1
%       n         norm    n-1      n-2
%------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef cheby_real < handle
    
    properties (Access = public)
        order           % order of polynomial expansion
        precision       % used for truncating the polynomial series
    end
    
    properties (Access = private)
        automatic       % toggle automatic determination of order
    end
    
    properties (Access = private)
        main_coeffs
        main_phase
        
        sub_alpha
        sub_coeffs
        sub_phase
    end
        
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = cheby_real (my_ord, my_prec)
            if isempty (my_ord)
                obj.order = 0;
            else
                obj.order = my_ord;
            end
            if isempty (my_prec)
                obj.precision = eps;
            else
                obj.precision = my_prec;
            end
            obj.automatic = false;
        end
        
        % Display propagator, overloading default disp method
        function disp(obj)
            global time
            prt.disp ('Chebychev polynomial expansion in real time')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp (['Kosloff number dE*dt / (2*hbar) : ' num2str(time.steps.m_alpha)])
            if obj.automatic
                prt.disp (['Automatic truncation if coefficients drop below : ' num2str(obj.precision)])
                prt.disp (['Number of polynomials required : ' int2str(obj.order)])
            else
                prt.disp (['Truncating polynomial series after : ' int2str(obj.order)])
            end
        end

        % Initialize propagator
        function init (obj, ~)
            global hamilt time
            
            % Set/check parameters
            if isfield(time,'pulse')
                prt.error ('This propagator is only for time independent Hamiltonians (no electric fields)')
            end
            
            if time.steps.m_alpha <= 40
                prt.error ('Kosloff number should not be below 40, then it becomes inefficient')
            end
     
            % Automatically determine expansion order
            if obj.order == 0
                obj.automatic = true;
                
                % Search backward
                for ii = round(2.0*time.steps.m_alpha) : -1 : 0.5*round(time.steps.m_alpha)
                    if abs(besselj(ii,time.steps.m_alpha))>obj.precision
                        break
                    end
                end
                obj.order = ii;
            end
            
            % Use Bessel functions to get expansion coefficients
            obj.main_coeffs = besselj(0:obj.order,time.steps.m_alpha);
            obj.main_coeffs(2:obj.order+1) = 2 * obj.main_coeffs(2:obj.order+1);
            
            % Get coefficients (columns) also for time sub steps (rows)
            obj.sub_alpha =  (1:time.steps.s_number)' * time.steps.s_delta * hamilt.range.delta/2;
            
            % Suddenly(?), BESSELJ  no longer accepts mixtures of row and column vectors
            [NU1,Z1] = meshgrid (0:obj.order,obj.sub_alpha);
            obj.sub_coeffs  = besselj(NU1,Z1);
            obj.sub_coeffs (:,2:obj.order+1) = 2 * obj.sub_coeffs(:,2:obj.order+1);
            
            % Phase factors to compensate for normalization of Hamiltonian
            obj.main_phase = exp( -1i * (hamilt.range.delta/2+hamilt.range.e_min) *                            time.steps.m_delta );
            obj.sub_phase  = exp( -1i * (hamilt.range.delta/2+hamilt.range.e_min) * (1:time.steps.s_number)' * time.steps.s_delta );
            
        end

        % Perform propagation
        function propa (obj, psi)
            
            global hamilt space time

            % Pre-allocate
            cheby0 = cell(hamilt.coupling.n_eqs,1);
            cheby1 = cell(hamilt.coupling.n_eqs,1);
            cheby2 = cell(hamilt.coupling.n_eqs,1);
            
            % Indexing of substeps, needed for autocorrelation
            n1 = time.steps.offset + 1;
            n2 = time.steps.offset + time.steps.s_number;
            
            %-----------------------------------------------------------
            %  Zero-th Chebyshev polynomial : phi_0 = 1
            %-----------------------------------------------------------
            for m = 1:hamilt.coupling.n_eqs
                cheby0{m} = psi.dvr{m};
                psi.sum{m} = obj.main_coeffs(1) * cheby0{m};
            end
            
            % Autocorrelation
            overlap = 0;
            for m = 1:hamilt.coupling.n_eqs
                overlap = overlap + sum ( conj(psi.ini{m}(:)).*cheby0{m}(:) .* space.weight(:) );
            end
            time.steps.acf(n1:n2) = obj.sub_coeffs(:,1) * overlap;
            
            %-----------------------------------------------------------
            %  First Chebychev polynomial phi_1 = -i*Hnorm
            %-----------------------------------------------------------
            apply_ham(psi,[0 0],1);
            for m = 1:hamilt.coupling.n_eqs
                cheby1{m} = - 1i  * psi.new{m};
                psi.sum{m} = psi.sum{m} + obj.main_coeffs(2) * cheby1{m};
            end
            
            % Autocorrelation
            overlap = 0;
            for m = 1:hamilt.coupling.n_eqs
                overlap = overlap + sum ( conj(psi.ini{m}(:)).*cheby1{m}(:) .* space.weight(:) );
            end
            time.steps.acf(n1:n2) = time.steps.acf(n1:n2) + obj.sub_coeffs(:,2) * overlap;
            
            %-----------------------------------------------
            %  Higher Chebychev polynomials (n>1) by recursion:
            %  phi_n = -2*i*Hnorm phi_{n-1} + phi_{n-2}
            %-----------------------------------------------
            for k=2:obj.order
                
                for m = 1:hamilt.coupling.n_eqs
                    psi.dvr{m} = cheby1{m};
                end
                apply_ham(psi,[0 0],1);
                for m = 1:hamilt.coupling.n_eqs
                    cheby2{m} = - 2 * 1i * psi.new{m} + cheby0{m};
                    psi.sum{m} = psi.sum{m} + obj.main_coeffs(k+1) * cheby2{m};
                    
                    cheby0{m} = cheby1{m};
                    cheby1{m} = cheby2{m};
                end
                
                % Autocorrelation
                overlap = 0;
                for m = 1:hamilt.coupling.n_eqs
                    overlap = overlap + sum ( conj(psi.ini{m}(:)).*cheby2{m}(:) .* space.weight(:) );
                end
                time.steps.acf(n1:n2) = time.steps.acf(n1:n2) + obj.sub_coeffs(:,k+1) * overlap;
                
            end
            
            %  Multiply wave function with complex phase factor
            for m = 1:hamilt.coupling.n_eqs
                psi.dvr{m} = psi.sum{m} * obj.main_phase;
            end
            
            %  Multiply autocorrelation with complex phase factor
            time.steps.acf(n1:n2) = time.steps.acf(n1:n2) .* obj.sub_phase;
            
        end
    end
end

