%------------------------------------------------------------------------------
%
% This class deals with coherent superpositions of complex-valued Gaussian
% wavepackets of given width, centered at r_0 / k_0 in position / momentum
% space. Wavefunction not yet normalized. We adopt a correlated description,
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

classdef gauss_corr < handle
    
    properties (Access = public)
        pos_0       % Matrix of the mean positions
        mom_0       % Matrix of the mean momenta
        width       % Matrix of the widths of the Gaussian packets
        factor      % Array of coefficients Cg for the Gaussian packets
    end
    
    properties (Access = private)
        n           % Number of Gaussians
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = gauss_corr
            
        end
        
        
        % Initialize Gaussian: Set/check parameters
        function init (obj)
            
            
            % Number of Gaussian wavepackets (number of rows in parameter matrix)
            obj.n = size ( obj.pos_0, 1 );
            
            % Set default values
            if isempty (obj.mom_0)
                obj.mom_0 = zeros(size(obj.pos_0));
            end
            if isempty (obj.factor)
                obj.factor = ones(obj.n, 1);
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
            prt.disp ( '                                                 ')
            prt.disp ( '           N     DOF      [                    ( Rk - R0gk )^2 ] ' )
            prt.disp ( 'psi(R) =  Sum C  Prd  exp [ I K0gk*(Rk-R0gk) - (-----------)   ] ' )
            prt.disp ( '          g=1  g k=1      [                    (   2*W_gk  )   ] ' )
            prt.disp ( '                                                 ')
            prt.disp ( 'Note that this wavepacket corresponds to a       ')
            prt.disp ( 'coherent (Glauber) state of a harmonic           ')
            prt.disp ( 'oscillator V(R) = k/2 (R-R0)^2 with mass m and   ')
            prt.disp ( 'width parameter W = (4*k*m)^(-1/4)               ')
            for g=1:obj.n
                prt.disp (' ')
                if obj.n>1
                    prt.disp ( [ 'Gaussian packet #' int2str(g)])
                    prt.disp (' ')
                    prt.disp ( [ 'Prefactor                 C    : ' num2str(     obj.factor(g))  ] )
                end
                prt.disp ( [ 'Mean value of position      R0 : ' num2str(     obj.pos_0(g,:)) ] )
                prt.disp ( [ 'Mean value of momentum      K0 : ' num2str(     obj.mom_0(g,:)) ] )
                prt.disp ( [ 'Position uncertainties       W : ' num2str(     obj.width(g,:)) ] )
                prt.disp ( [ 'Momentum uncertainties 1/(2*W) : ' num2str(0.5./obj.width(g,:)) ] )
            end
        end
        
        
        
        % Evaluate wave function on a grid
        function wave (obj, psi)
            
            global space
            
            psi.dvr{1} = zeros(size(space.dvr{1}));
            
            % Loop over Gaussian packets and sum them up
            for g=1:obj.n
               
                gauss = ones(size(space.dvr{1}));
                
                for k = 1:space.n_dim
                    gauss = gauss .* exp (  1i*(space.dvr{k}-obj.pos_0(g,k)) *  obj.mom_0(g,k)...
                        -((space.dvr{k}-obj.pos_0(g,k)) / (2*obj.width(g,k))).^2  );
                end

                psi.dvr{1} = psi.dvr{1} + obj.factor(g) * gauss;
                
            end
        end
    end
end

