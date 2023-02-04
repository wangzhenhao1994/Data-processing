% The Legendre DVR is used for problems where the minor quantum number (most often
% called m) is conserved while the major quantum number (typically named l) can
% vary. Note that for homonuclear diatomic molecules, symmetry arguments require
% that l is always either even or odd. If your potential does not fulfill this
% automatically (linearly polarized lasers do), you have to add the required
% symmetry. Note that the Legendre DVR considers the quantum number l as momentum.
%        
% Rotational quantum numbers from l = m_0 to (and including) l_max, 
% therefore it has (l_max - m_0 + 1) elements which is also the 
% number of grid points in position space.


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%               2017 Burkhard Schmidt
%
% see the README file for license details.

classdef legendre < handle
    
    properties (Access = public)
        
        label       % labelling the degree of freedom
        dof         % which degree of freedom
        mass        % mass thet enters the kinetic energy 
        nokin       % enable/disable the kinetic energy operator
 
        R_dof       % Index of the degree of freedom that is used for the R
        R_0         % Fixed value of R for rigid rotor

                
        weight      % weights in DVR ("position space")
        x_grid      % grid points in DVR ("position space")
        p_grid      % grid points in FBR ("momentum space")

        l_max       % maximum angular momentum of the spectral basis
        m_0         % (constant) value of the minor quantum number. 

        
    end
    
    properties (Access = private)
        
        kin         % DVR grid representation of the kinetic energy operator
        kin_expo    % same as grid, but exponentiated, used with the split operator 
        
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
        function obj = legendre
            obj.dof = 1;
            obj.mass = 1;
            obj.m_0 = 0;
            obj.nokin = false;
        end
        
        function disp (obj) % Diplay method, overloading default disp method
            % Works only if init_grid has been called previously
            prt.disp ( '***************************************************************' )
            prt.disp ( [ 'DVR for the degree of freedom: ' obj.label] )
            prt.disp ( 'Discretization scheme: Gauss-Legendre DVR' )
            prt.disp ( '***************************************************************' )
            prt.disp ( ' ' )
            prt.disp ( 'Parameters of the (angular) grid' )
            prt.disp ( [ 'Number of grid points       : ' int2str(obj.l_max - abs(obj.m_0) + 1) ] )
            prt.disp ( [ 'Lower bound for L           : ' int2str(abs(obj.m_0)) ] )
            prt.disp ( [ 'Upper bound for L           : ' int2str(obj.l_max) ] )
            prt.disp ( [ 'Azimuthal quantum number m  : ' int2str(obj.m_0) ] )
            prt.disp ( ' ' )
        end
        
        % Evaluate kinetic energy function (for trajectories!)
        function eval_kin (obj, state)
            prt.error ('Code missing')
        end
        
    end
    
    methods (Access = private) % in separate files within *private* directory
        
        [shapedmat, permutation, shapedims] = shape(obj, input)
                 
    end

    methods (Access = private, Static = true)

        [weights, xi, trafo] = quadrature (num_points, m)
        reshaped_mat = shape_back(input, permutation, shapedims)

    end
    
end
