%--------------------------------------------------------------------------
%
% Generic superclass from which all 
% other main classes can inherit:
% wave, traj, ket, rho, 
% sht, mssh, fssh, sssh, ...
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-.... Burkhard Schmidt
%
% see the README file for license details.

classdef generic < handle
    
    properties (Access = public)
        
        string0     % Name of this class (or rather, its subclasses)
        
        save_export % Toggle saving to (binary) files
        save_dir    % Directory where to save|load
        save_file   % File name template
        save_suffix % File name suffix
        save_step   % Start a new file after this many steps
        save_mem    % Max. file size: 500 MB
        file_name   % Fiule name (with suffix and extension)
        
    end
    
    methods (Access = public)
        
        % Constructor: Default values for file I/O
        function obj = generic

            obj.save_export = false;     % Toggle saving to (binary) files
            obj.save_dir    = pwd;       % Directory where to save|load 
            obj.save_file   = [];        % File name template
            obj.save_suffix = [];        % File name suffix
            obj.save_step   = [];        % Start a new file after this many steps 
            obj.save_mem    = 500*2^20;  % Max. file size: 500 MB
   
        end
        
        % More methods: see separate files
        save_0 ( obj )                   % Save general setup
        save_1 ( obj )                   % Save (first time step only)
        load_0 ( obj )                   % Save general setup

    end
    
end

