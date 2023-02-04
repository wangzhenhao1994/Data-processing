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
%               2008-2010 Ulf Lorenz
%
% see the README file for license details.

function save_0 ( state )

% Setting default file name template (=class name) unless specified otherwise
if isempty(state.save_file)
    state.save_file = state.string0;
end

% Contruct file name from template and (optional) suffix
if isempty (state.save_suffix)
    state.file_name = strcat(state.save_file, '_0.mat');
else
    state.file_name = strcat(state.save_file, '_', state.save_suffix, '_0.mat');
end

% Console & logfile output
prt.disp('***************************************************************');
prt.disp('Saving general information to mat-file:')
prt.disp('Expectation values, Hamiltonian, representations etc')
prt.disp('***************************************************************');
prt.disp(' ');
prt.disp(['Directory          : ' prt.shorten(state.save_dir) ])
prt.disp(['File name          : ' prt.shorten(state.file_name) ])
prt.disp(['Steps per mat-file : ' int2str(state.save_step) ])
prt.disp(' ');

% Make directory if necessary
if ~isfolder(state.save_dir)
    mkdir(state.save_dir);
end

