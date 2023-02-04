% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef metiu2 < pot.generic & handle
    
    properties (Access = public)
        
        d_e         % dissociation energy
        r_e         % equilibrium position
        alf         % range parameter
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = metiu2
        end

        % Initialize potential: Set/check parameters
        function init (~)
            global hamilt space
            if space.n_dim ~= 2
                prt.error('This potential works only for two DOFs')
            end
            
            if ~isa(space.dof{1}, 'dof.fft') || ~isa(space.dof{2}, 'dof.legendre')
                prt.error('This potential works only for one FFT and one Legendre DOF')
            end
            
            if hamilt.coupling.n_eqs > 1
                prt.error('This potential works only for one Schroedinger equation')
            end
        end
    
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            prt.disp('Model potential for angular problem ')
            prt.disp('***************************************************************')
            prt.disp('                                                  ')
            prt.disp(' Dateo & Metiu, J.Chem.Phys. 95, 7392 (1991)      ')
            prt.disp('                                                  ')
            prt.disp('                                                 2')
            prt.disp(' V(r,\Theta) = V (r) + 0.1*D*(1-cos \Theta)*(r /r)')
            prt.disp('                M                             e   ')
            prt.disp('                                                  ')
            prt.disp([ 'Dissociation energy d_e:' num2str(obj.d_e)])
            prt.disp([ 'Range parameter a      :' num2str(obj.alf)])
            prt.disp([ 'Equilibrium dist. r_e  :' num2str(obj.r_e)])
            prt.disp (' ')
        end

        % Evaluate grid representation of potential energy functions
        function V = V(obj,r)
            
            % Morse potential
            V = obj.d_e * (1 - exp(-obj.alf ...
                * (r{1} - obj.r_e))).^2;
            
            % Add the cosine \theta term
            V = V + obj.d_e * 0.1 ...
                * (1 - r{2}) .* (obj.r_e ./ r{1}).^2;
            
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            prt.error ('Code for calculation of forces still missing')
        end
        
    end
end
