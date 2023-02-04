%------------------------------------------------------------------------------
%
% Wigner (phase space) plots of 1-dim wavepacket dynamics:
% Grid representations of total energy and Wigner functions
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%
% see the README file for license details.

function wigner ( state, step )
global expect hamilt plots space

% Wigner transform for FFT grids only
if ~isa (space.dof{1}, 'dof.fft')
    prt.error ('Cannot make a Wigner plot of a non-Fourier grid')
end

% Grid representation of total energy function(s)
if plots.density.energy && step==1
    
    % Create grid representation of total energy function compatible
    % Use only inner part of kinetic energy grid and double size
    kinetic = space.dof{1}.kin;
    kinetic = kinetic(space.dof{1}.n_pts/4+1 : 3*space.dof{1}.n_pts/4);
    kinetic = [kinetic kinetic]';
    kinetic = reshape ( kinetic, space.dof{1}.n_pts, 1 );

    for m = 1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.pot{m,m}.dvr)
            potential = hamilt.pot{m,m}.dvr;
        else
            potential = zeros(size(space.dvr{1}));
        end
        hamilt.tef.wig {m}  = ...
            repmat ( potential, 1, space.dof{1}.n_pts)' + ...
            repmat (   kinetic, 1, space.dof{1}.n_pts);

        hamilt.tef.wig{m}( hamilt.tef.wig{m}>hamilt.truncate.e_max ) = hamilt.truncate.e_max;
        hamilt.tef.wig{m}( hamilt.tef.wig{m}<hamilt.truncate.e_min ) = hamilt.truncate.e_min;

    end
end

if isa(state,'traj') % this includes sub-classes for SHT
    return
end

% Perform Wigner transform of all wavefunctions
state.wig_max = 0;
for m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step) > expect.min_pop
        state.wig{m} = wave.wigner ( state.dvr{m}' );
        state.wig_max = max ( state.wig_max, max(abs(state.wig{m}(:))   ) );
    end
end
