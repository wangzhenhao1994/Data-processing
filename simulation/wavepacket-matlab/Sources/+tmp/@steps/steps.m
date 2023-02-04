%------------------------------------------------------------------------------
%
% Temporal discretization: Sub-steps
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2010-2011 Ulf Lorenz
%
% see the README file for license details.

classdef steps < handle
    
    properties (Access=public)
        
        m_start     % index of first main time step
        m_stop      % index of last main time step
        m_delta     % size of main time steps
        m_number    % number of main time steps
        m_grid      % vector containing the main time steps 
        m_alpha     % Kosloff number for main time steps

        t_total     % total simulation time
        
        s_number    % number of sub steps per main time step
        s_delta     % size of substeps
        s_grid      % vector containing the substeps
        s_alpha     % Kosloff number for substeps
        
        offset      % Offset of substeps in main propagation step
        
        acf         % autocorrelation function
        freq        % frequencies
        spec        % spectrum
        
    end
    
    methods (Access=public)
        
        function obj = steps
            
            obj.m_start  = 0;            % by default start main steps from zero
            obj.m_stop   = 1;            % by default one propagation step
            obj.m_delta  = 1;            % default size of main time steps
            obj.s_number = 1;            % by default no sub-stepping
            
        end
        
        function init (obj)
            
            global hamilt
            
            % Main steps
            obj.m_number = abs(obj.m_stop - obj.m_start) + 1;              % number of time steps
            obj.t_total  = obj.m_delta * ( obj.m_number - 1 );                  % total simulation time
            obj.m_grid   = (obj.m_start + [ 0 : obj.m_number-1 ]') * obj.m_delta; % Column(!) vector
            
            % Sub steps
            obj.s_delta  = obj.m_delta / obj.s_number;       
            obj.s_grid   = obj.m_start*obj.m_delta ...
                + [0:(obj.m_number-1)*obj.s_number]' * obj.s_delta;
            
            prt.disp ('***************************************************************')
            prt.disp ('Temporal discretization: Main steps and sub steps')
            prt.disp ('***************************************************************')
            prt.disp ( ' ' )
            prt.disp ( [ 'Index of first step  : ' int2str(obj.m_start)   ] )
            prt.disp ( [ 'Index of last step   : ' int2str(obj.m_stop)    ] )
            prt.disp ( [ 'Number of steps      : ' int2str(obj.m_number-1)] )
            prt.disp ( ' ' )
            prt.disp ( [ 'Size of steps        : ' num2str(obj.m_delta)   ] )
            prt.disp ( [ 'Total time           : ' num2str(obj.t_total)   ] )
            prt.disp ( ' ' )
            prt.disp ( [ 'Number of sub-steps  : ' int2str(obj.s_number)  ] )
            prt.disp ( [ 'Size of sub-steps    : ' num2str(obj.s_delta)   ] )
            prt.disp ( ' ' )
            
            % hamilt.range.delta will be only set/used when this function
            % is called from qm_propa but not from qm_optimal
            if isfield (hamilt,'range')
                if isfield (hamilt.range,'delta')
                    obj.m_alpha = obj.m_delta * hamilt.range.delta / 2;
                    obj.s_alpha = obj.s_delta * hamilt.range.delta / 2;
                    prt.disp ( [ 'Spectral range of Hamiltonian : ' num2str(hamilt.range.delta)  ] )
                    prt.disp ( [ 'Kosloff alpha (main steps)    : ' num2str(obj.m_alpha)  ] )
                    prt.disp ( [ 'Kosloff alpha (sub steps)     : ' num2str(obj.s_alpha)  ] )
                    prt.disp ( ' ' )
                end
            end
            
        end
        
        spectrum ( obj )
            
    end
        
end
