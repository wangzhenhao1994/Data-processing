%**********************************************************************************
% Kinetic operator in FFT degree of freedom acting on wavefunction 
% ----------------------------------------------------------------
%
% To apply the operator, the following three steps are performed:
% (1) The wavefunction is expanded in the spectral basis 
%     here: IFFT transform to obtain FBR
% (2) The kinetic operator is applied in this representation
% (3) The wavefunction is converted back to pseudospectral representation
%     here: FFT transform to obtain DVR
%
% Original wavefunction is assumed in field psi.dvr or psi.new
% (if new is true)
% Resulting wavefunction is stored in psi.new
%**********************************************************************************

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function kinetic ( obj, psi, new )

global hamilt

% If grid_kinetic is called without previous initialisation, initialise it here.
if isempty(obj.kin)
	init_kin(obj, 1);
end

if obj.nokin
    for m = 1:hamilt.coupling.n_eqs
        psi.new{m} = zeros(size(psi.dvr{m}));
    end
    return
end

% Note: This is the same as grid_dvr2fbr + multiplication with obj.kin +
% grid_fbr2dvr. Without normalisation, fftshift, and excess function calls,
% it should be somewhat faster.
for m = 1:hamilt.coupling.n_eqs
    if new
        psi.new{m} = ifft( obj.kin_shift .* fft(psi.new{m}, [], obj.dof), [], obj.dof);
    else
        psi.new{m} = ifft( obj.kin_shift .* fft(psi.dvr{m}, [], obj.dof), [], obj.dof);
    end
end
