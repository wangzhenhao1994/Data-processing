% ------------------------------------------------------------------------
%
% Loads the setup saved in previous simulations from an external file.
% The path and filename of saved calculation to be loaded are taken from
% properties of object "state" provided as the first parameter in calling.
%
% If a second parameter "toggle" is supplied (boolean with true value),
% this function will also overwrite the following properties of class
% "state" with those from the saved calculation:
% expect, hamilt, info, space, time
% Global variables that are not overwritten: plot.
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20xy Burkhard Schmidt
%
% see the README file for license details.

function load_0 ( state, choice )

global expect hamilt info space time

% Inherit from method of superclass
load_0@generic ( state );

% Optionally set global variables
if choice > 0
    
    load(fullfile(state.save_dir, state.file_name), ...
        'M_h', 'M_a', 'M_d', 'M_p', 'M_s', ...
        'M_m', 'M_v', 'M_l', 'M_o', '-mat');

    state.M_ham = M_h;
    state.M_amo = M_a;
    state.M_dip = M_d;
    state.M_pol = M_p;
    state.M_sbc = M_s;
    state.M_mat = M_m;
    state.M_vec = M_v;
    state.M_lab = M_l;
    state.M_obs = M_o;
    
    prt.disp('Setting matrix representations');

end

if choice > 1
    
    load(fullfile(state.save_dir, state.file_name), ...
        'expect', 'hamilt', 'info', 'space', 'time', '-mat')
    
    prt.disp('Setting global variables: expect, hamilt, info, space, time');
    
end