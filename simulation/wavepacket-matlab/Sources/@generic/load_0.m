% ------------------------------------------------------------------------
%
% Loads the setup saved in previous simulations from an external file.
% The path and filename of saved calculation to be loaded are taken from
% properties of object "state", which can be of class "traj" (trajectories)
% or of class "wave" (wavefunctions) or any of the other main classes and
% which is provided as an input argument.
%
% The class property "save_step" indicates how many time steps are stored
% together in each of the mat-files. This information is needed, e.g. in
% "load" which actually retrieves the trajectories or wavefunctions or etc.
% (class "traj" or "wave" or etc.) from those data files.
%
% If a second parameter "toggle" is supplied (boolean with true value),
% this function will also overwrite the following global variables with
% those from the saved calculation: expect, hamilt, info, space, time
% Global variables that are not overwritten: plot.
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20xy Burkhard Schmidt
%               2008-2010 Ulf Lorenz
%
% see the README file for license details.

function load_0 ( state )

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

% Get variable "save_step" from data file
load(fullfile(state.save_dir, state.file_name), 's_s', '-mat');
state.save_step = s_s;

% Console & logfile output
prt.disp('***************************************************************');
prt.disp('Loading general information from mat-file:')
prt.disp('Expectation values, Hamiltonian, representations etc')
prt.disp('***************************************************************');
prt.disp(' ');
prt.disp(['Directory name     : ' prt.shorten(state.save_dir)  ])
prt.disp(['File name template : ' prt.shorten(state.file_name) ])
prt.disp(['Steps per mat-file : ' int2str(state.save_step) ])
prt.disp(' ');

end