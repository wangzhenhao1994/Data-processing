%--------------------------------------------------------------------------
% Potential energy: Paired quantum dot (1D) for one or two electrons 
% Note that both singlet and triplets are included
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2020 Federico Pont
%
% see the README file for license details.

classdef pqd1D < pot.generic & handle
    
    properties (Access = public)
        
        b_L         % width of left QD (inverse square of ...)
        V_L         % depth of left QD
        b_R         % width of right QD (inverse square of ...)
        V_R         % depth of right QD
        R           % distance between QDs
        lambda      % Coulomb interaction
                
    end
    
    methods (Access = public)
        
       % Constructor: Set default values
        function obj = pqd1D
            obj.V_L = 0;
            obj.V_R = 0;
            obj.lambda = 0;
        end
        
        
        % Initialize potential: Set/check parameters
        function init (obj)
            global space
            
            if (space.n_dim > 2)
                prt.error ('Should only be used in 1D for one or two electron simulations')
            end
            
            if (obj.b_L<=0) || (obj.b_R<=0) || (obj.R<=0) 
                prt.error ('The quantum dots sizes and distance should be positive')
            end
            
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            prt.disp ('Paired quantum dot 1D potential for one or two electrons')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp ( [ 'Inverse square of left dot size    : ' num2str(obj.b_L) ] )
            prt.disp ( [ 'Left dot depth                     : ' num2str(obj.V_L) ] )
            prt.disp ( [ 'Inverse square of Right dot size   : ' num2str(obj.b_R) ] )
            prt.disp ( [ 'Right dot depth                    : ' num2str(obj.V_R) ] )
            prt.disp ( [ 'Distance between QDs               : ' num2str(obj.R)   ] )
            prt.disp (' ')
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)
            
            % Summing up contributions from each component of position vector
            V = zeros(size(r{1}));
            
            % Loop over (one or two) electrons
            for k = 1:length(r)
                V = V + obj.V_L ...
                    * ( exp ( -obj.b_L * ( r{k}+obj.R/2 ).^2 ) ) ...
                    + obj.V_R ...
                    * ( exp ( -obj.b_R * ( r{k}-obj.R/2 ).^2 ) );
            end
            
            % Coulomb repulsion
            if length(r)>1 && obj.lambda ~= 0 
                V = V + obj.lambda*sqrt(pi/2)*erfcx ( abs(r{1}-r{2})/sqrt(2) );
            end
            
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(~,~)
            prt.error ('Code missing')
        end
    end
end

