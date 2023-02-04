%--------------------------------------------------------------------------
%
% Compose figure with two or three subplots
% displaying various expectation values as time series
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

classdef expect < vis.styles & handle
    
    properties (Access = public)
        
        w_left     % Left border of plot window
        w_lower    % Left border of plot window
        w_width    % Width of plot window
        w_height   % Height of plot window

        errorbar    % Toggle plotting errorbars
        legends     % Toggle curve legends
        logo        % Toggle logos in all four corners
        wide        % Toggle wide format: 16:9
        
        p_min       % Lower bound for population
        p_max       % Upper bound for populations
        
        e_min       % Lower bound for energies
        e_max       % Upper bound for energies
        
        export      % Toggle graphics export
        file        % Set the output to a custom filename (w/o suffix)
        
    end
    
    properties (Access = private)
        
        mask_tot    % Fake: Find steps already existing
        mask_cha    % Find steps where populations of individual channels are non-negligible

    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = expect
            
            % Inheriting from superclass
            obj = obj@vis.styles;
            
            % Position and size of plot window: Format 7:9
            obj.w_left   = round(11*obj.s_height/16);  % Left border
            obj.w_lower  = round(01*obj.s_height/16);  % Lower border
            obj.w_width  = round(07*obj.s_height/16);  % Width
            obj.w_height = round(09*obj.s_height/16);  % Height
            obj.w_width =  2*round(obj.w_width /2);    % Should be even
            obj.w_height = 2*round(obj.w_height/2);    % Should be even
            
            % Appearence of plot
            obj.errorbar = false;        % Toggle plotting errorbars
            obj.legends  = true;         % Toggle legends
            obj.logo     = true;         % Toggle logos in all four corners
            obj.wide = [];               % Toggle wide format: 16:9
            
            % Range settings for populations
            obj.p_min = -0.1;            % Minimum of population plot
            obj.p_max = +1.1;            % Maximum of population plot
            
            % Range settings for energies
            obj.e_min   = [];            % Lower bound of the energy plot
            obj.e_max   = [];            % Upper bound of the energy plot
            
            % Export to graphics file: jpg/tiff/epsc/pdf, etc ...
            obj.export  = false;         % Toggle graphics export
            obj.file    = [];            % Set custom filename (suffix determines image type)
            
            % Masks where individual/total populations do exist
            obj.mask_tot = [];           % Fake: Find steps already existing
            obj.mask_cha = [];           % Find steps where individual populations 
           
        end
        
        %--------------------------------------------
        % Show expectation values
        %--------------------------------------------
        function show_expect (obj,~)
            
            global expect hamilt info state time
            
            if ~isa(state,'wave') && ~isa(state,'traj')
                prt.error('"Expect" plots available only for wavefunctions or trajectory bundles')
            end
          
            if obj.wide % Wide format: 16:9
                w=16; h=09; o=9/16;
            else % Narrow format: 7:9
                w=07; h=09; o=0;
            end
            
            % Find steps where populations of individual channels are non-negligible
            for m=1:hamilt.coupling.n_eqs
                obj.mask_cha{m} = find ( expect.pop.cha{m} > expect.min_pop );
            end
            
            % Fake: Find steps already existing
            obj.mask_tot = find ( expect.pop.tot > expect.min_pop );
            
            % Time dependent simulations
            switch lower (info.program)
                case {'qm_propa'}
                    
                    if isfield(time,'pulse')
                        % External electric field versus time
                        subplot ( 'Position',[o+1/w 1/h 5/w 7/(3*h)] )
                        cla;
                        efield (obj)
                        if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    elseif ~isempty(time.steps.acf)
                        % Autocorrelation function versus time
                        subplot ( 'Position', [o+1/w 1/h 5/w 7/(3*h)] )
                        cla;
                        correlation (obj)
                        if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    end
                
                    % Expectation values: Energies versus time
                    subplot ( 'Position', [o+1/w (1+7/3)/h 5/w 7/(3*h)] )
                    cla;
                    energies (obj)
                    if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    
                    % Populations/projections vs. time
                    subplot ( 'Position', [o+1/w (1+14/3)/h 5/w 7/(3*h)] )
                    cla;
                    population (obj)
                    if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    
                    % Time independent simulations
                case {'qm_bound'}
                    
                    % Populations/projections and energies
                    if isfield(hamilt, 'amo') || hamilt.coupling.n_eqs>1
                        subplot ( 'Position', [o+1/w 1/h 5/w 3/h] )
                        energies (obj)
                        if obj.wide; set ( gca, 'YAxisLocation', 'right'); end

                        subplot ( 'Position', [o+1/w 5/h 5/w 3/h] )
                        population (obj)
                        if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    else
                        % Energies only
                        subplot ( 'Position', [o+1/w 1/h 5/w 7/h] )
                        energies (obj)
                        if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    end
                    
                otherwise
                    prt.error ('Wrong choice of program')
                
            end
            
        end
 
 
        % Details of surface hopping trajectory methodology
        function show_detail (obj,state)
            global info time
            
            if obj.wide % Wide format: 16:9
                w=16; h=09; o=9/16;
            else % Narrow format: 7:9
                w=07; h=09; o=0;
            end
            
            if isa (state,'wave')
                return
            end
            
            pos_siz = [o+1/w 1.5/h 5/w 1/h];  % x y w h   
           
            h = annotation('textbox');
            set(h, ...
                'Position',            pos_siz, ...
                'HorizontalAlignment', 'left', ...
                'Interpreter',         'none', ...
                'LineStyle',           'none', ...
                'String',              {state.string1, state.string2, state.string3, state.string4, state.string5, state.string6, time.propa.string}, ...
                'FontName',            obj.f_name, ...
                'FontSize',            obj.f_large, ...
                'FontWeight',          obj.f_heavy)
            if strcmpi(info.system,'Matlab')
                set (h,'FitHeightToText', 'on')
            elseif strcmpi(info.system,'Octave')
                set(h,'FitBoxToText', 'on')
            end
                        
        end
        
    end
    
    
    methods (Access = private)
        
        population  (obj)
        energies     (obj)
        correlation (obj)
        efield      (obj)
        
    end
    
end
