% ------------------------------------------------------------------------
%
% Stores "ket" (state vector) data in internal arrays and 
% finally saves them to MATLAB® formatted file 
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019      B. Schmidt
%               2008-2010 Ulf Lorenz
%
% see the README file for license details.

function save ( state, step )

global time

persistent arr_x index dim

% If we do not want to save the trajectories, exit immediately
if ~state.save_export
    return;
end
    

% First time step only
if step == 1

    % Length of state vectors
    dim = length(state.x);
    
    % Initialize saving
    save_1 (state)  % from superclass "generic"
    
    % Create arrays for states to be saved
    arr_x = zeros(dim, state.save_step);
    
    % Initialize counting of states
    index = 1;
    
end

% Save the ket vectors to array and increase index
arr_x(:,index) = state.x;
index = index + 1;

% Save the arrays in a Matlab formatted file
if index > state.save_step || step == time.steps.m_number

    % Get file name (containing index)
	file_indx = int32(floor( (step-1)/double(state.save_step) ))+1; 
    if isempty (state.save_suffix)
        file_name = strcat(state.save_file, '_', int2str(file_indx), '.mat');
    else
        file_name = strcat(state.save_file, '_', state.save_suffix, '_', int2str(file_indx), '.mat');
    end
    prt.disp('***************************************************************');
    prt.disp('Saving state vectors to mat-file')
    prt.disp('***************************************************************');
    prt.disp(' ');
    prt.disp(['Directory : ' prt.shorten(state.save_dir) ])
    prt.disp(['File name : ' prt.shorten(file_name) ])
    prt.disp(' ');
    
    % Save arrays to data file and reset them afterwards
    save(fullfile(state.save_dir,file_name), 'arr_x', '-mat');
    arr_x = zeros(dim, state.save_step);
    index = 1;
    
end

% Last time step only
if step == time.steps.m_number
    save_0 (state)
end

end
