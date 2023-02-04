%--------------------------------------------------------------------------
%
% Close coupling scheme and diabatic/adiabatic 
% representation of the Hamiltonian operator
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef coupling < handle
    
    properties (Access=public)
        
        n_eqs       % Number of coupled equations
        labels      % Labels of coupled equations
        represent   % Diabatic or adiabatic
        
        ini_rep     % Adiabatic or diabatic initial state
        ini_coeffs  % Coefficients of initial state (adiabatic or diabatic)
        ini_norm    % Toggle normalization of initial state

    end
    
    methods (Access=public)
        
        % Constructor: Set default values
        function obj = coupling
            obj.n_eqs = 1;
            obj.labels = [];
            obj.represent = [];
            obj.ini_norm = true;
            obj.ini_rep = 'dia';
        end
        
        % Initialization
        function init (obj, state)
            
            % Labels of (coupled) states
            if isempty(obj.labels)
                for m = 1:obj.n_eqs
                    obj.labels{m} = int2str(m);
                end
            end
            
            % Default: diabatic representation
            if isempty (obj.represent)
                obj.represent='dia';
            end
            
            % For single channel wavefunctions: always adiabatic
            if isa(state,'wave') && obj.n_eqs == 1
                obj.represent='adi';
            end
            
            % Normalize vector of coefficients
            if ~isempty(obj.ini_coeffs) && obj.ini_norm
                obj.ini_coeffs = obj.ini_coeffs / norm (obj.ini_coeffs);
            end
       
        end
        
        % Display properties of object
        function disp (obj)
            
            if obj.n_eqs == 1
                return
            end
            
            prt.disp ('***************************************************************')
            prt.disp ('Initialize close coupling scheme    ')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            switch lower (obj.represent)
                case 'dia'
                    prt.disp ( 'Diabatic representation' )
                case 'adi'
                    prt.disp ( 'Adiabatic representation' )
                otherwise
                    prt.error ( [ 'Incorrect choice of representation  : ' obj.represent ] )
            end
            
            prt.disp ( [ 'Number of (coupled) equations      : ' int2str(obj.n_eqs) ] )
            for m=1:obj.n_eqs
                prt.disp ( [ int2str(m) ': ' obj.labels{m} ] )
            end
            prt.disp (' ')
            
            prt.disp ('***************************************************************')
            prt.disp ('Initial distribution of populations')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp ( [ '(A)diabatic representation    : '         obj.ini_rep        ] )
            prt.disp ( [ 'Normalize initial state       : ' int2str(obj.ini_norm     ) ] )
            prt.disp (' ')
            
            if  ~isempty(obj.ini_coeffs)
                prt.disp ( [ 'Coefficients of wavefunctions : ' num2str(obj.ini_coeffs   ) ] )
                prt.disp ( [ 'Populations of wavefunctions  : ' num2str(abs(obj.ini_coeffs).^2) ] )
                
                % Check validity of coefficient vector
                if any (obj.ini_coeffs>1)
                    prt.error('Initial coefficients must not exceed one')
                end
            end
            
            prt.disp (' ')
            
        end
        
    end
    
end

