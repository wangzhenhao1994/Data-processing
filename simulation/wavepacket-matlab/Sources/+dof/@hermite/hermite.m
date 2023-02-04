% The Hermite DVR is suited for an expansion in eigenfunctions of the harmonic
% oscillator (Hermite polynomials times an exponential). In contrast to other DVR
% methods, the grid points are not natively bounded to some interval, but can lie
% anywhere (of course, they are always located somewhere around the minimum of the
% harmonic oscillator).
% 
% Furthermore, the Gauss-Hermite quadrature is conveniently defined in scaled
% coordinates (corresponding to m * omega = 1, r_e = 0).  To use them for real
% problems, you have to supply the shift and the properties of the harmonic
% potential you want to use.

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%               2017 Burkhard Schmidt
%
% see the README file for license details.

classdef hermite < handle
    
    properties (Access = public)
        
        label       % labelling the degree of freedom
        dof         % which degree of freedom
        mass        % mass that enters the kinetic energy
        nokin       % enable/disable the kinetic energy operator
        
        omega       % angular frequency of the potential: k = m * omega^2.
        v_2         % force constant of Harmonic oscillator properties
        
        n_pts       % number of grid points (== number of basis functions)
        r_e         % equilibrium position of the harmonic oscillator.
        
        weight      % weights in DVR ("position space")
        x_grid      % grid points in DVR ("position space")
        p_grid      % grid points in FBR ("momentum space")
        
    end
    
    properties (Access = private)
        
        kin         % DVR grid representation of the kinetic energy operator
        kin_expo    % Same as grid, but exponentiated for split operator
        kin_max     % Maximum of kinetic energy
        dvr_mom     % DVR(!) grid representation of the momentum operator
        
        trafo2fbr   % transformation matrix : DVR=>FBR
        trafo2dvr   % transformation matrix : FBR=>DVR
        
    end
    
    methods (Access = public) % in separate files within *this* directory
        
        fbr = dvr2fbr(obj, dvr)
        dvr = fbr2dvr(obj, fbr)
        init_grid(obj)
        init_kin(obj, fraction, output)
        kinetic     (obj, psi, new)
        kinetic_exp (obj, psi)
        dvrkin = kinetic2dvr(obj)
        dvr = matrix2dvr(obj, fbr)
        retval = momentum(obj, psi)
        obj = subsasgn(obj, index, val)
        retval = subsref(obj, index)
        
        % Constructor method: Initialize default property values
        function obj = hermite
            obj.dof   = 1;
            obj.mass  = 1;
            obj.omega = 1;
            obj.r_e   = 0;
            obj.nokin = false;
            obj.kin_max = 0;
        end
        
        function disp (obj) % Diplay method, overloading default disp method
            % Works only if init_grid has been called previously
            prt.disp ( '***************************************************************' )
            prt.disp ( [ 'DVR for the degree of freedom: ' obj.label] )
            prt.disp ( 'Discretization scheme: Gauss-Hermite DVR' )
            prt.disp ( '***************************************************************' )
            prt.disp ( ' ' )
            prt.disp ( 'Parameters of the grid' )
            prt.disp ( ['Number of grid points : ' int2str(obj.n_pts)])
            prt.disp ( ['scaling factor        : ' num2str(sqrt(obj.mass*obj.omega))])
            prt.disp ( ['harmonic frequency    : ' num2str(obj.omega)])
            prt.disp ( ['position of minimum   : ' num2str(obj.r_e)  ])
            prt.disp ( ' ' )
        end
        
        % Evaluate kinetic energy function (for trajectories!)
        function eval_kin (obj, state)
            prt.error ('Code missing')
        end
        
    end

    methods (Access = private)
        [shapedmat, permutation, shapedims] = shape(obj, input)
    end
    
    methods (Access = private, Static = true)
        
        [weights, xi] = quadrature (num_points)
        reshaped_mat = shape_back(input, permutation, shapedims)
                
    end
    
end
