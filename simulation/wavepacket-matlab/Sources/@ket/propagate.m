%--------------------------------------------------------------------------
%
% Propagate ket object (state vector) subject to a given Hamiltonian
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
    time.steps.offset = + 1;
    init ( time.propa, obj );
    prt.disp ('***************************************************************')
    prt.disp ('Numeric propagation scheme:                       ')
    disp ( time.propa )
    prt.disp (' ')
else
    % Actually do the propagation, using adaptive sub-stepping
    time.steps.offset = (step-2) * time.steps.s_number + 1;
    propa ( time.propa, obj, step );   
end
