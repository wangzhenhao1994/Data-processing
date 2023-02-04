% ------------------------------------------------------------------------
%
% Loads the setup saved in previous simulations from an external file.
% The path and filename of saved calculation to be loaded are taken from
% properties of object "state" provided as the first parameter in calling.
%
% If a second parameter "toggle" is supplied (boolean with true value),
% this function will also overwrite the following properties of class
% "state" with those from the saved calculation: ABNCD etc.
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20xy Burkhard Schmidt
%
% see the README file for license details.

function load_0 ( state, choice )

global info time

% Inherit from method of superclass
load_0@generic ( state );

% Optionally set global variables
if choice > 0
    
    load(fullfile(state.save_dir, state.file_name), ...
        'A', 'B', 'N', 'C', 'D', ...
        'Q', 'S', 'T', 't', ...
        'x_i', 'x_e', 'x_l', ...
        'y_i', 'y_e', 'y_l', '-mat');
    
    state.A = A;
    state.B = B;
    state.N = N;
    state.C = C;
    state.D = D;
    state.Q = Q;
    state.S = S;
    state.T = T;
    state.title = t;
    state.x_initial = x_i;
    state.x_equilib = x_e;
    state.x_label   = x_l;
    state.y_initial = y_i;
    state.y_equilib = y_e;
    state.y_label   = y_l;
    
    prt.disp('Setting properties of state: ABNCD etc and x_i, x_e, x_l, y...');
    
end

if choice > 1
    
    load(fullfile(state.save_dir, state.file_name), ...
        'info', 'time', '-mat');
    
    prt.disp('Setting global variables: info, time');
    
end


end