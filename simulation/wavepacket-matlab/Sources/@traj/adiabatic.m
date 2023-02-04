%--------------------------------------------------------------------------
%
% Transform (matrix-valued) potential given at the grid points
% from diabatic to adiabatic representation.
%
% Note that within the context of classical and quantum-classical WavePacket 
% simulations, grid-based representations of the potential energy are used 
% only for visualization, hence only for space.n_dim<3
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

function adiabatic (~,step,~)
global hamilt space

% Use this method only *before* propagation
if step >= 0
    return
end

% Use this method only in 1D or 2D
if space.n_dim>2
    return
end

% Use this method only if adiabatic representation is desired
if strcmpi(hamilt.coupling.represent,'dia')
    return
end

% Get adiabatic potentials by solving eigenproblem
switch hamilt.coupling.n_eqs
    case 1 % Single channel
        ham.pot_eig_1 ( 0, 0 )
        
    case 2 % Two channels
        ham.pot_eig_2 ( 0, 0 )
        
    otherwise % Multiple channels
        ham.pot_eig_N ( 0, 0 )       
end

% Save adiabatic potentials
for m = 1:hamilt.coupling.n_eqs
    hamilt.pot{m,m}.dvr = hamilt.eig_val{m}; % diagonal only
end

end