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

% Do not save potentially huge arrays from qm_bound
if strcmp(info.program,'qm_bound')
    hamilt.eigen.matrix    = [];
    hamilt.eigen.eig_vecs  = [];
    hamilt.eigen.transform = [];
end

% Apparently, Matlab cannot save class properties
s_s = state.save_step;
M_h = state.M_ham;
M_a = state.M_amo;
M_d = state.M_dip;
M_p = state.M_pol;
M_s = state.M_sbc;
M_m = state.M_mat;
M_v = state.M_vec;
M_l = state.M_lab;
M_o = state.M_obs;

save(fullfile(state.save_dir, state.file_name), 's_s', ...
    'expect', 'hamilt', ...
    'info', 'plots', ...
    'space', 'time', ...
    'M_h', 'M_a', 'M_d', 'M_p', 'M_s', ...
    'M_m', 'M_v', 'M_l', 'M_o', '-mat');
end
