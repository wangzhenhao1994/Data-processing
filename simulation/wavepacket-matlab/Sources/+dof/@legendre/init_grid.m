% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007 Martin Winter
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function init_grid(obj)

% Gaussian quadrature
[obj.weight, obj.x_grid, obj.trafo2fbr] = obj.quadrature(obj.l_max - abs(obj.m_0) + 1, abs(obj.m_0));

% The spectral grid is sqrt( l(l+1) )
obj.p_grid = (abs(obj.m_0):obj.l_max)';
obj.p_grid = sqrt(obj.p_grid.*(obj.p_grid+1));

% The inverse trafo is just the transpose since it is unitary (even orthogonal)
obj.trafo2dvr = obj.trafo2fbr';

% Add the weights to the transformation matrix, since it
% involves an integration/summation over the DVR grid.
obj.trafo2fbr = obj.trafo2fbr .* repmat(obj.weight', [ length(obj.weight) 1] );
