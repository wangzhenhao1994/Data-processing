% ------------------------------------------------------------------------
%
% Loads states (wavefunctions!)
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

global hamilt

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
prt.disp('Loading wave functions from mat-file')
prt.disp('***************************************************************');
prt.disp(' ');
prt.disp(['Directory/folder   : ' prt.shorten(state.save_dir)  ])
prt.disp(['File name template : ' prt.shorten(file_name) ])
prt.disp(['Index of this step : ' int2str(step_indx) ])
prt.disp(' ');

% Load wavefunction arrays from data file
load(fullfile(state.save_dir,file_name), 'arr_dvr', '-mat');
for m=1:hamilt.coupling.n_eqs
    state.dvr{m} = arr_dvr{m,step_indx};
end

end
