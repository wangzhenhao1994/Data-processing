%--------------------------------------------------------------------------
%
% Propagate trajectories subject to a given Hamiltonian
% using one of the user-specified ODE solvers.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008,2010 Ulf Lorenz
%
% see the README file for license details.

function propagate (obj, step)
global time

% Initialize one of the propagator classes of your choice
if step==1
    eq_init (time.eq_motion, obj, time.propa)
    save_previous(obj)
else
    % Loop over substeps: Call propagator classes of your choice
    for k = 1:time.steps.s_number
        eq_propa (time.eq_motion, obj, time.propa); 
        save_previous(obj);
    end
    
    % Get new energies
    apply_ham ( obj );
    
end
