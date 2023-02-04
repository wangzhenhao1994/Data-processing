% ------------------------------------------------------------------------
%
% Saves additional information (Hamiltonian etc) in a
% MATLAB formatted binary file (MAT-file, ....mat) 
% Should be called after the *last* time step only
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019      B. Schmidt
%
% see the README file for license details.

function save_0 ( state )

global expect hamilt info plots space time

% Inherit from method of superclass
save_0@generic ( state );

% Apparently, Matlab cannot save class properties
s_s = state.save_step;

save(fullfile(state.save_dir, state.file_name), 's_s', ...
    'expect', 'hamilt', ...
    'info', 'plots', ...
    'space', 'time',  '-mat');

end
