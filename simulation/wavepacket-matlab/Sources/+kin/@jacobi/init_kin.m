% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008,2011 Ulf Lorenz
%
% see the README file for license details.

function init_kin (obj, fraction, output)

% initialises the kinetic energy

global hamilt space time

if nargin < 3
    output = true;
end

%% Checks

if isempty( obj.mass_R ) || isempty( obj.mass_r )
     prt.error ('masses not set')
end
if isempty( obj.dof_c ) 
     prt.error ('c degree of freedom not set')
end
if isempty( obj.dof_R ) 
     prt.error ('R degree of freedom not set')
end
if isempty( obj.dof_r ) && isempty( obj.r_0 )
     prt.error ('r degree of freedom not set')
end

% Check that we have a reasonable grid (i.e. Legendre polynomials)
if ~isa ( space.dof{obj.dof_c}, 'dof.legendre' )
     prt.error ('Jacobi kinetic energy only works for Legendre grids')
end


%% Informational output

if output
    prt.disp ('***************************************************************')
    prt.disp ('Kinetic energy operator: Triatomic molecule ABC                ')
    prt.disp ('***************************************************************')
    prt.disp ('                                                  ')
    prt.disp ('           [     1             1    ] ^ 2         ')
    prt.disp (' T (c) = - [ --------  +   -------- ] L           ')
    prt.disp ('           [ 2 M R^2       2 m r^2  ]             ')
    prt.disp ('                                                  ')
    prt.disp ('where R, M is the distance and reduced mass of BC,')
    prt.disp ('where r is the distance of A to BC center of mass,')
    prt.disp ('where m is the reduced mass of A and BC.          ')
    prt.disp ('See, e.g.,  J. Chem. Phys. 116, 4403 (2002)       ')
    prt.disp ('                                                  ')
    prt.disp ( [ 'M     : ' num2str(obj.mass_R) ] )
    prt.disp ( [ 'm     : ' num2str(obj.mass_r) ] )
    prt.disp ( [ 'DOF c : ' num2str(obj.dof_c)  ] )
    prt.disp ( [ 'DOF R : ' num2str(obj.dof_R)  ] )
    if ~isempty(obj.dof_r)
        prt.disp ( [ 'DOF r : ' num2str(obj.dof_r) ] )
    else
        prt.disp ( [ 'constant r : ' num2str(obj.r_0) ] )
    end
    prt.disp (' ')
end


%% Create all the grids

% Prefactor
if isempty( obj.r_0 )
    obj.grid = 1 ./ (2 * space.dvr{obj.dof_R}.^2 * obj.mass_R) ...
               + 1 ./ (2 * space.dvr{obj.dof_r}.^2 * obj.mass_r);
else
    obj.grid = 1 ./ (2 * space.dvr{obj.dof_R}.^2 * obj.mass_R) ...
               + 1 ./ (2 * obj.r_0.^2 * obj.mass_r);
end

% L^2
obj.grid = space.fbr{obj.dof_c}.^2 .* obj.grid;

% Truncation
obj.grid(obj.grid > hamilt.truncate.delta) = hamilt.truncate.delta;

% Short-time propagator
obj.grid_exp = exp(-1i * obj.grid * time.steps.s_delta * fraction);
