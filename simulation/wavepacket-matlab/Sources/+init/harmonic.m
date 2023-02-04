%------------------------------------------------------------------------------
%
% This class creates the initial state from a (ground
% or excited) eigenstate of a harmonic oscillator.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2009 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%
% see the README file for license details.

classdef harmonic < init.generic & handle
    
    properties (Access = public)

        m_r         % reduced mass
        r_e         % equilibrium distance
        n_q         % quantum number

        omega       % harmonic frequency
        v_2         % force constant
                
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = harmonic
            obj.r_e = 0;
            obj.n_q = 0;
        end
        
        % Initialize pendular functions: Set/check parameters
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
                prt.error ('Quantum number must not be negative')
            end
            
            % Use either v_2 or omega from the input file; calculate missing quantity
            if isfield(obj, 'v_2') && ~isfield(obj, 'omega')
                obj.omega = sqrt(obj.v_2 / obj.m_r);
            end
            if isfield(obj, 'omega') && ~isfield(obj, 'v_2')
                obj.v_2 = obj.omega^2 * obj.m_r;
            end
            if isfield(obj, 'omega') && isfield(obj, 'v_2')
                prt.error ('Specify either omega or v_2 but not both of them')
            end
            
        end

       % Display Gaussian, overloading default disp method
       function disp(obj)
           prt.disp ('Harmonic oscillator eigenstate')
           prt.disp ('***************************************************************')
           prt.disp ('   ' )
           prt.disp ('Harmonic oscillator eigenstate')
           prt.disp ( ['Angular Frequency            : ' num2str(obj.omega)] )
           prt.disp ( ['corresp. Force constant      : ' num2str(obj.v_2)] )
           prt.disp ( ['(Reduced) mass               : ' num2str(obj.m_r)] )
           prt.disp ( ['Equilibrium position         : ' num2str(obj.r_e)] )
           prt.disp ( ['Quantum number (eigenstate)  : ' int2str(obj.n_q)] )
       end
           
        % Evaluate wave function on a grid ==> Q/M propagation
        function wave (obj)
            global space

            % HO eigenstate \Psi(x) = exp(-m\omega/2 * x^2) * H_n(\sqrt(m\omega) * x)
            factor    = obj.m_r * obj.omega;
            position  = space.dvr{obj.dof} - obj.r_e;
            obj.dvr = exp(- factor/2 * position.^2) ...
                .* math.hermite(sqrt(factor) * position, obj.n_q);
            
        end
        
        % Sample phase space density ==> Q/C propagation
        function traj (obj, traj)
            prt.error ('Code for phase space sampling still missing')
        end
        
    end
    
end
