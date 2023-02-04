%------------------------------------------------------------------------------
%
% Stationary states of a simple planar pendulum
%
%      V(R) = V0/2 [ 1 + cos (m*(R-R0) ] 
%
% This class gives the eigenfunctions of a simple planar pendulum 
% which can be given in terms of ordinary (angular) Mathieu function. 
% To be used as initial wavefunction (defined on a ring) for quantum pendula 
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%
% see the README file for license details.

classdef pendulum1 < init.generic & handle
    
    properties (Access = public)
        
        multiple    % multiplicity; should be 1 or 2
        order       % order of elliptic function
        barrier     % height of the barriers
        shift       % angular shift
        parity      % parity: 'c' or 's' oir 'l'

    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = pendulum1
            
        end
        
        % Initialize pendular functions: Set/check parameters
        function init (obj)
                        
            switch obj.multiple
                case 1 % Only even integer orders r
                    if obj.order~=2*round(obj.order/2)
                        prt.error ('Order of Mathieu functions not compatible with m!')
                    end
                case 2 % Only integer orders r
                    if obj.order~=round(obj.order)
                        prt.error ('Order of Mathieu functions not compatible with m!')
                    end
                otherwise
                    prt.error ('Fractional order Mathieu functions not yet implemented!')
            end

        end
        
        % Display pendular functions, overloading default disp method
        function disp(obj)
            
            prt.disp ('Pendular state wavefunction')
            prt.disp ('***************************************************************')
            prt.disp ('                                                  ')
            prt.disp ('Ordinary (angular) Mathieu functions              ')
            prt.disp ('                                                  ')
            prt.disp ('ce_r (eta,q) with r>=0 : cosine elliptic functions')
            prt.disp ('se_r (eta,q) with r>=1 :   sine elliptic functions')
            prt.disp ('                                                  ')
            prt.disp ('solving a Schrödinger equation for potential      ')
            prt.disp ('                                                  ')
            prt.disp ('       V(R) = V0/2 [ 1 + cos (m*(R-R0) ]          ')
            prt.disp ('                                                  ')
            prt.disp ('for periodic coordinate 0 < R <= 2*pi where the   ')
            prt.disp ('multiplicity m indicates number of maxima/minima  ')
            prt.disp ('and where eta=m*(R-R0)/2 and q=2V0/m^2            ')
            prt.disp ('                                                  ')
            prt.disp ('m=1: Only pi-periodic solutions: Only even r      ')
            prt.disp ('m=2: All 2pi-periodic solutions: Integer r        ')
            prt.disp ('m>2: Requires also fractional orders of r         ')
            prt.disp ('     Yields complex solutions: not yet implemented')
            prt.disp ('  ')
            
            switch lower(obj.parity)
                case ('c')
                    prt.disp ( 'Cosine elliptic' )
                 case {'s'}
                    prt.disp ( 'Sine elliptic' )
                 case {'l'}
                    prt.disp ( 'Localized combination of cosine and sine elliptic' )
                otherwise
                    prt.error ('Wrong choice for parity of Mathieu function')
            end
            prt.disp ( [ 'Mathieu order r      : ' int2str( obj.order)    ] )
            prt.disp ( [ 'Multiplicity m       : ' int2str( obj.multiple) ] )
            prt.disp ( [ 'Potential barrier V0 : ' num2str( obj.barrier)  ] )
            prt.disp ( [ 'Angular shift R0     : ' num2str( obj.shift)    ] )
            
        end
        
        % Evaluate wave function on a grid  ==> Q/M propagation
        function wave (obj)
            global space
            
            % Arguments of Mathieu equation
            r = obj.order;
            q = obj.barrier*2/obj.multiple^2;
            e = obj.multiple*(space.dvr{obj.dof}-obj.shift)/2;
            
            switch lower(obj.parity)
                case ('c')
                    obj.dvr = math.cee(r,q,e);
                case {'s'}
                    obj.dvr = math.see(r,q,e);
                case {'l'}
                    obj.dvr = (math.cee(r,q,e)+math.see(r+1,q,e))/sqrt(2);

            end           
        end
        
        % Sample phase space density ==> Q/C propagation
        function traj (obj, traj)
            prt.error ('Code for phase space sampling still missing')
        end        
        
    end
end
            
            
            
