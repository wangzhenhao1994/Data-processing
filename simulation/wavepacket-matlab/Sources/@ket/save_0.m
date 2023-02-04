% ------------------------------------------------------------------------
%
% Saves additional information (ABNCD etc) in a
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

global info plots time

% Inherit from method of superclass
save_0@generic ( state );

% Apparently, Matlab cannot save class properties
s_s = state.save_step;
A = state.A;
B = state.B;
N = state.N;
C = state.C;
D = state.D;
Q = state.Q;
S = state.S;
T = state.T;
t = state.title;
x_i = state.x_initial;
x_e = state.x_equilib;
x_l = state.x_label;
y_i = state.y_initial;
y_e = state.y_equilib;
y_l = state.y_label;

save(fullfile(state.save_dir, state.file_name), 's_s', ...
    'info', 'plots', 'time', ...
    'A', 'B', 'N', 'C', 'D', ...
    'Q', 'S', 'T', 't', ...
    'x_i', 'x_e', 'x_l', ...
    'y_i', 'y_e', 'y_l', '-mat')

end
