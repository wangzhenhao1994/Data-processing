% ------------------------------------------------------------------------
%
% Loads states (trajectories!) 
% from an external mat-file.
%
% To be done: Cache the above arrays, i.e. make them persistent.
% Upon every subsequent call check if we already cached this file, 
% otherwise load and cache it
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008-2010 Ulf Lorenz
%
% see the README file for license details.

function load ( state, step )

global space

% Calculate indices for the file suffix and the position in the file.
file_indx = int32(floor( (step-1)/double(state.save_step) ))+1;
step_indx = step - (file_indx-1)*state.save_step;

% Get file name (containing index)
if isempty (state.save_suffix)
    file_name = strcat(state.save_file, '_', int2str(file_indx), '.mat');
else
    file_name = strcat(state.save_file, '_', state.save_suffix, '_', int2str(file_indx), '.mat');
end
prt.disp('***************************************************************');
prt.disp('Loading trajectories from mat-file')
prt.disp('***************************************************************');
prt.disp(' ');
prt.disp(['Directory/folder   : ' prt.shorten(state.save_dir)  ])
prt.disp(['File name template : ' prt.shorten(file_name) ])
prt.disp(['Index of this step : ' int2str(step_indx) ])
prt.disp(' ');

% Load trajectory arrays from data file
load(fullfile(state.save_dir,file_name), 'arr_pos', 'arr_mom', 'arr_cha', '-mat');
state.cha = arr_cha{1,step_indx};
for k=1:space.n_dim
    state.pos{k} = arr_pos{k,step_indx};
    state.mom{k} = arr_mom{k,step_indx};
end

end
