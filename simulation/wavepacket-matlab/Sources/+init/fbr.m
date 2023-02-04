%------------------------------------------------------------------------------
%
% This function creates an initial state as a "momentum", or rather FBR eigenstate,
% (single state or coherent superposition of states)
% where the exact nature of the state depends on the grid used in the DVR scheme.
% For FFT grids, it returns a plane wave, for Legendre polynomials, the result
% is a spherical harmonic, for Hermite grids it is a Harmonic oscillator eigenstate.
%
% However, since you have to specify that you want the n-th eigenfunction, this
% routine is somewhat awkward to use for plane waves, as their quantum numbers
% translate into momenta going from -pmax to +pmax. I.e. to get momentum 0, you
% have to specify something like "#grid points / 2".
%
% You can either set up a pure eigenstate by specifying the variable state or a
% superposition by specifying the array of coefficients in the variable coeffs.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% see the README file for license details.

classdef fbr  < init.generic & handle

    properties (Access = public)
        state       % index of a single(!) FBR state
        coeffs      % Coefficients for several FBR states

    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = fbr
        end

        % Initialize fbr function: Set/check parameters
        function init (obj)
            if ~isempty(obj.state) && ~isempty(obj.coeffs)
                prt.error ('Specify either "state" or "coeffs" but not both!')
            end
           if isempty(obj.state) && isempty(obj.coeffs)
                prt.error ('Either "state" or "coeffs" have to be specified!')
            end
        end

        % Display fbr function, overloading default disp method
        function disp(obj)
            
            prt.disp ('FBR ("momentum") grid eigenfunction.               ')
            prt.disp ('***************************************************************')
            prt.disp (' ' )
            if ~isempty(obj.state) 
                prt.disp ( ['Selecting the ' int2str(obj.state) '-th eigenstate'])

            end
            if ~isempty(obj.coeffs)
                prt.disp( ['Using coefficient vector : ' num2str(obj.coeffs) ])
            end
            
        end

        % Evaluate wave function on a grid ==> Q/M propagation
        function wave (obj)
            global space
            
            % Single state
            if ~isempty(obj.state)    
                obj.dvr = zeros(size(space.fbr{obj.dof}));
                help1 = obj.state;
                help2 = space.dof{obj.dof}.p_grid;
                help3 = help2(help1);
                obj.dvr( space.fbr{obj.dof} == help3 ) = 1;
            end
                
            % Coherent superposition 
            if ~isempty(obj.coeffs)
                obj.dvr = zeros(size(space.fbr{obj.dof}));
                for ii = 1:length(obj.coeffs)
                    help2 = space.dof{obj.dof}.p_grid;
                    help3 = help2(ii);
                    obj.dvr(space.fbr{obj.dof} == help3 ) = obj.coeffs(ii);
                end
            end
            
            % Transform the result to DVR space
            obj.dvr = fbr2dvr(space.dof{obj.dof}, obj.dvr);

        end

        % Sample phase space density ==> Q/C propagation
        function traj (obj, traj)
            prt.error ('Code for phase space sampling still missing')
        end
        
    end
end
