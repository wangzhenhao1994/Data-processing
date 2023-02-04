%------------------------------------------------------------------------------
%
% This class deals with coherent superpositions of complex-valued Gaussian
% wavepackets of given width, centered at r_0 / k_0 in position / momentum
% space. Wavefunctions not yet normalized, for one degree of freedom only.
%
% The function parameter tells us which degree of freedom we are supposed to
% fill up. The resulting wave function is returned by this function.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% see the README file for license details.

classdef gauss < init.generic & handle

    properties (Access = public)
        pos_0       % Vector of the mean positions
        mom_0       % Vector of the mean momenta
        width       % Vector of the widths of the Gaussian packets
        factor      % Vector of coefficients Cg for the Gaussian packets
    end
    
    properties (Access = private)
        n           % Number of Gaussians
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = gauss
            
        end
        
        
        % Initialize Gaussian: Set/check parameters
        function init (obj)
            
            % Number of Gaussian wavepackets
            obj.n = size ( obj.pos_0, 1 );
            
            % Set default values
            if isempty ( obj.mom_0 )
                obj.mom_0 = zeros(obj.n,1);
            end
            if isempty (obj.factor)
                obj.factor = ones(obj.n,1);
            end
            
        end
        
        % Display Gaussian, overloading default disp method
        function disp(obj)         
            if obj.n==1
                prt.disp ( 'Gaussian minimum uncertainty wave packet' )
            else
                prt.disp ( ['Coherent superposition of ' int2str(obj.n) ' Gaussian packets'])
            end
            prt.disp ('***************************************************************')
            prt.disp ( '                                                       ' )
            prt.disp ( '           N         [                 ( R - R0g )^2 ] ' )
            prt.disp ( 'psi(R) =  Sum C  exp [ I K0g*(R-R0g) - (---------)   ] ' )
            prt.disp ( '          g=1  g     [                 (  2*Wg   )   ] ' )
            prt.disp ( '                                                       ' )
            prt.disp ( 'Note that this wavepacket corresponds to a coherent    ' )
            prt.disp ( '(Glauber) state of a harmonic oscillator with potential' )
            prt.disp ( 'V(R) = k/2 (R-R0)^2 = 1/2 m omega^2 (R-R0)^2 with mass ' )
            prt.disp ( 'm, frequency omega = (k/m)^(1/2), and width            ' )
            prt.disp ( 'W = (4*k*m)^(-1/4) = (2*m*omega)^(-1/2)                ' )
            for g=1:obj.n
                prt.disp (' ')
                if obj.n>1
                    prt.disp ( [ 'Gaussian packet #' int2str(g)])
                    prt.disp (' ')
                    prt.disp ( [ 'Prefactor                 C  : ' num2str(     obj.factor(g))] )
                end
                prt.disp ( [ 'Mean value of position    R0 : ' num2str(     obj.pos_0(g)) ] )
                prt.disp ( [ 'Mean value of momentum    K0 : ' num2str(     obj.mom_0(g)) ] )
                prt.disp ( [ 'Position uncertainty       W : ' num2str(     obj.width(g)) ] )
                prt.disp ( [ 'Momentum uncertainty 1/(2*W) : ' num2str(0.5./obj.width(g)) ] )
            end
         end
        
        % Evaluate wave function on a grid ==> Quantum propagation
        function wave (obj)

            global space
            
            obj.dvr = zeros(size(space.dvr{obj.dof}));
            
            % Loop over Gaussian packets and sum them up
            for g=1:obj.n
                
                obj.dvr = obj.dvr + obj.factor(g) * ...
                    exp (  1i*(space.dvr{obj.dof}-obj.pos_0(g)) *  obj.mom_0(g)...
                    -((space.dvr{obj.dof}-obj.pos_0(g)) / (2*obj.width(g))).^2  );
            end
        end
        
        % Sample phase space density ==> Q/C propagation
        function traj (obj, traj)
            
            % Choose Gaussian packet randomly
            g = randi (obj.n,traj.n_p, 1);
            
            % Normally distributed random positions
            if traj.n_p==1
                traj.pos {obj.dof} = obj.pos_0(g);
            else
                traj.pos {obj.dof} = obj.width(g) .* randn([traj.n_p 1]) + obj.pos_0(g);
            end
            
            % Normally distributed random momenta
            if traj.n_p==1
                traj.mom {obj.dof} = obj.mom_0(g);
            else
                traj.mom {obj.dof} = 0.5./obj.width(g) .* randn([traj.n_p 1]) + obj.mom_0(g);
            end
            
        end
        
    end
end
