%--------------------------------------------------------------------------
%
% Verlet back-propagation
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-.... Leonardo Cancissu Araujo
%
% see the README file for license details.

function [pos, mom, frc, pot_mats, frc_mats, D_1, U_1] = verlet_prop_back(obj,pos,mom,frc,n)
global time space

% Propagate momenta by half time step
for d = 1:space.n_dim
    mom{d} = mom{d} - 0.5 * time.steps.s_delta .* frc{d};
end

% Propagate positions by full time step
for d = 1:space.n_dim
    pos{d} = pos{d} - time.steps.s_delta .* mom{d} / space.dof{d}.mass;
end

% Calculate and save all diabatic potentials in advance
pot_mats = ham.pot_dia(pos);

% Compute adiabatic potential matrix and eigenvector matrix
[U_1 , D_1] = ham.pot_eig_adi(pot_mats);

% Calculate and save all diabatic forces in advance
frc_mats = ham.frc_dia(pos);

frc = ham.frc_adi ( pot_mats, frc_mats , U_1 , n );

% Propagate momenta by half time step
for d = 1:space.n_dim
    mom{d} = mom{d} - 0.5 * time.steps.s_delta .* frc{d};
end
end
