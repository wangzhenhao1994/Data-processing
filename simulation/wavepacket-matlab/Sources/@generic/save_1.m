% ------------------------------------------------------------------------
%
% Prepare for saving states (traj, wave, etc) in an  
% MATLAB® formatted binary file (MAT-file, ....mat) 
% Should be called after the *first* time step only
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019      B. Schmidt
%               2008-2010 Ulf Lorenz
%
% see the README file for license details.

function save_1 ( state )

global time

% Setting default file name template (=class name) unless specified otherwise
if isempty(state.save_file)
    state.save_file = state.string0;
end

% Output of file/directory
prt.disp('***************************************************************');
prt.disp('Saving states to mat-files')
prt.disp('***************************************************************');
prt.disp(' ');
prt.disp(['Directory name     : ' prt.shorten(state.save_dir)    ])
prt.disp(['File name template : ' prt.shorten(state.save_file)   ])
prt.disp(['File name suffix   : ' prt.shorten(state.save_suffix) ])
if ~isempty(state.save_step)
    prt.disp(['Maximum steps per file    : ' int2str(state.save_step) ])
end
if state.save_mem>2^30
    prt.disp(['Maximum file size  : ' num2str(state.save_mem/2^30) ' GB' ])
elseif state.save_mem>2^20
    prt.disp(['Maximum file size  : ' num2str(state.save_mem/2^20) ' MB' ])
elseif state.save_mem>2^10
    prt.disp(['Maximum file size  : ' num2str(state.save_mem/2^10) ' kB' ])
else
    prt.disp(['Maximum file size  : ' int2str(state.save_mem     ) ' B' ])
end

% Memory per time step:
mem_per_step = memory_size (state);
if mem_per_step>2^30
    prt.disp(['Memory per step    : ' int2str(mem_per_step/2^30) ' GB' ])
elseif mem_per_step>2^20
    prt.disp(['Memory per step    : ' int2str(mem_per_step/2^20) ' MB' ])
elseif mem_per_step>2^10
    prt.disp(['Memory per step    : ' int2str(mem_per_step/2^10) ' kB' ])
else
    prt.disp(['Memory per step    : ' int2str(mem_per_step     ) ' B' ])
end

% Determine the stepsize;
steps_per_file = min(int32(ceil(state.save_mem / mem_per_step)), ...
    time.steps.m_number);

% Optionally: User defined step size
if isempty(state.save_step) || state.save_step > steps_per_file
    state.save_step = steps_per_file;
end
prt.disp(['Steps per mat-file : ' int2str(state.save_step) ])
prt.disp(' ');

