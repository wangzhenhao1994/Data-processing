%------------------------------------------------------------------------
%
% This class deals with bound eigenstates of a Morse oscillator.
% Not using the standard recursion to generate Laguerre polynomials as im-
% plemented in math.laguerre but using instead the modified version from 
% J. P. Dahl, M. Springborg,  J. Chem. Phys. 88, 4535-47 (1988)
%
% Potential function
%
%      V(r) = d_e*(1-exp(-alf*(r-r_e)))^2
%
% m_r     - (reduced) mass
% d_e     - dissociation energy
% r_e     - equilibrium position
% alf     - range parameter (steepness)
% n_q     - quantum number of the desired eigenfunction
%
%------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2009 Burkhard Schmidt
%               2009 Ulf Lorenz
%
% see the README file for license details.

classdef morse < init.generic & handle
    
    properties (Access = public)
        
        d_e         % dissociation energy
        r_e         % equilibrium position
        alf         % range parameter
        
        m_r         % reduced mass
        n_q         % quantum number
        
    end
    
    properties (Access = private)
        omg         % Harmonic frequency 
        chi         % Anharmonicity
        lbd         % lambda = omega / (2*chi)
        n_b         % Number of bound states
        e_n         % Energy of n_q(th) state
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = morse
            obj.r_e = 0;
            obj.n_q = 0;
        end
        
        % Initialize Morse functions: Set/check parameters
        function init (obj)
            global space
                   
           % Default: Mass from the grid setting along the respective degree of freedom 
           if ~isfield(obj, 'm_r')
               obj.m_r = space.dof{obj.dof}.mass;
           end

           % Check quantum number
           if round(obj.n_q)~=obj.n_q
               prt.error ('Quantum number should be integer')
           end
           if (obj.n_q<0)
                prt.error ('quantum number must not be negative')
            end
           
           % Check parameters
           if ( obj.d_e<=0 || obj.r_e<=0 || obj.alf<=0 || obj.m_r<=0 )
                prt.error ('Morse prameters should be positive')
           end
                      
           % Harmonic frequency and anharmonicity constant
           obj.omg = obj.alf * sqrt ( 2*obj.d_e/obj.m_r );
           obj.chi = obj.omg^2 / ( 4*obj.d_e );
           obj.lbd = obj.omg / (2*obj.chi); % = sqrt(2*d_e*m_r)/a_e
           
           % Number of bound states (counting from zero!)
           obj.n_b = floor (obj.lbd - 1/2);
           
           if (obj.n_q>obj.n_b)
               prt.error ('quantum number exceeds number of bound states')
           end
           
           % Bound state energy
           obj.e_n = ...
               + obj.omg*(obj.n_q+1/2) ...
               - obj.chi*(obj.n_q+1/2)^2 ...
               - obj.d_e;
           
        end
        
        % Display Gaussian, overloading default disp method
        function disp(obj)
            
            prt.disp ('Morse oscillator eigenstate')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp (['Dissociation energy     : ' num2str(obj.d_e)])
            prt.disp (['Equilibrium position    : ' num2str(obj.r_e)])
            prt.disp (['Range parameter (alfa)  : ' num2str(obj.alf)])
            prt.disp (['(Reduced) mass          : ' num2str(obj.m_r)])
            prt.disp (' ')
            prt.disp (['Number of bound states  : ' int2str(obj.n_b+1)])
            prt.disp (['Harmonic frequency      : ' num2str(obj.omg)])
            prt.disp (['Anharmonicity           : ' num2str(obj.chi)])
            prt.disp (' ')
            prt.disp (['Quantum number          : ' int2str(obj.n_q)])
            prt.disp (['Energy                  : ' num2str(obj.e_n)])
            
        end
        
        % Evaluate  wave function on a grid  ==> Q/M propagation
        function wave (obj)
            global space
 
            % Transformed arguments; parameters
            y  = obj.alf * ( space.dvr{obj.dof}-obj.r_e ); % Eq. 33
            xi = 2*obj.lbd * exp(-y);                           % Eq. 39
            s  = 2*obj.lbd - 2*obj.n_q - 1;                     % Eq. 46

            % Recursion
            psi.n = ones(size(space.dvr{obj.dof}));
            if (obj.n_q>0)
                psi.n1=zeros(size(psi.n));
                for n=1:obj.n_q
                    psi.n2=psi.n1;
                    psi.n1=psi.n;
                    psi.n = ((2*n+s-1-xi).*psi.n1 - sqrt((n-1)*(n+s-1))*psi.n2); % Eq. 49
                    psi.n = psi.n / sqrt(n*(s+n));
                end
            end

            % Normalization from Dahl's paper works only for nq=0 ?!?
            % psi.n = psi.n .* xi.^(s/2) .* exp(-xi/2) .* sqrt(alf*s*gamma(n_q+1)/gamma(2*lambda-n_q));

            % Tricky normalization (thanks to Marius Lewerenz, Marne le Vallee)
            scalog = ( log(xi)*s - xi + log(obj.alf*s) - gammaln(s+1) ) / 2 ;
            obj.dvr = psi.n .* exp(scalog);
            
        end
        
        % Sample phase space density ==> Q/C propagation
        function traj (obj, traj)
            prt.error ('Code for phase space sampling still missing')
        end        
        
    end
end