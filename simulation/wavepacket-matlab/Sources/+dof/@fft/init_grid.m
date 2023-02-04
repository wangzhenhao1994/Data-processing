% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function init_grid(obj)

% this is the (constant!) weight in coordinate space 
obj.x_dlt = (obj.x_max - obj.x_min) / obj.n_pts;
obj.weight = ones(obj.n_pts, 1) * obj.x_dlt;

% the grid points in coordinate space as a column (!) vector
obj.x_grid = linspace(obj.x_min, obj.x_max-obj.x_dlt, obj.n_pts)';

% The mapping from grid points in the DVR pseudospectral basis (coordinate space)
% to the FBR spectral basis (momentum space)
obj.p_min = - pi / obj.x_dlt;
obj.p_max = - obj.p_min;
obj.p_dlt = 2*obj.p_max / obj.n_pts;

% the grid points in momentum space as a column (!) vector
obj.p_grid = linspace(obj.p_min, obj.p_max - 2*obj.p_max/obj.n_pts, obj.n_pts)';