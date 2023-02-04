% Copyright (C) 2010 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space time params


% This file only works with the associated 'rundemo.m' script.  It uses a
% global variable "params" from which the time stepping, wave function
% initialisation, and wave function saving behaviour are copied. The rundemo.m
% script then recycles this initialise.m function to send the initial and target
% wavepackets back and forth in time.

% Number of coupled equations
hamilt.coupling.n_eqs  = 2;
hamilt.coupling.represent = 'dia';
hamilt.coupling.labels = {'X^1 A_1', 'A^1 B_1'};

% Grid definition; masses are not needed here since we use
% a custom operator; it would be smarter to use an unbalanced
% grid (O-D is heavier and thus needs more grid points than O-H),
% but we shall stick to the original setup as close as possible.
space.dof{1}       = dof.fft;            % using FFT grid
space.dof{1}.label = 'O-D';              % label for axis descriptions
space.dof{1}.n_pts = 512;                % number of points
space.dof{1}.x_min = 1.15;               % lower boundary of the grid
space.dof{1}.x_max = 26.75;              % upper boundary of the grid

space.dof{2}       = dof.fft;            % using FFT grid
space.dof{2}.label = 'O-H';              % label for axis descriptions
space.dof{2}.n_pts = 512;                % number of points
space.dof{2}.x_min = 1.15;               % lower boundary of the grid
space.dof{2}.x_max = 26.75;              % upper boundary of the grid

hamilt.amo{1} = amo.reaction;            % projecting on D + OH channel
hamilt.amo{1}.reac  = 1;
hamilt.amo{1}.prod  = 2;
hamilt.amo{1}.side  = 'r';
hamilt.amo{1}.label = 'D + OH';

hamilt.amo{2} = amo.reaction;            % projecting on H + OD channel
hamilt.amo{2}.reac  = 1;
hamilt.amo{2}.prod  = 2;
hamilt.amo{2}.side  = 'p';
hamilt.amo{2}.label = 'H + OD';

% Temporal discretisation copied from the params global variable.
% The electric field is also copied.
time = params.time;

% Hamiltonian operator
hamilt.truncate.e_min = -0.5;            % lower truncation of energy
hamilt.truncate.e_max = 0.5;             % upper truncation of energy

hamilt.kin{1}       = kin.triatomic;     % Kinetic energy for fixed bending angle
hamilt.kin{1}.theta = 104.52 * pi/180;   % bending angle
hamilt.kin{1}.dof   = [1 2];             % indices of coordinates for AB, BC distance
hamilt.kin{1}.mass  = [2.13 16.00 1.00] / atomic.m.u;
                                         % note that the D mass is off by a few percent
                                         % due to approximations in Ashwanis setup
                                         % (he assumed mu_OD = 2 * mu_OH).

hamilt.pot{1,1} = pot.H2O;               % potential energy (ground state)
hamilt.pot{2,2} = pot.H2O;               % potential energy (excited state)

hamilt.dip{1}      = cell(2);
hamilt.dip{1}{1,2} = dip.H2O;            % transition dipole moments

% Initial wave functionand some details on whether or not we save the wave
% function are also copied from the externally setup params variable.
psi = params.psi;

% Typically, we do not want to have plots, but let the rundemo script decide.
% If it does not prohibit plotting, there are a few standard settings we still
% can do.
if params.plot
    
    % Plot time evolution of density
    plots.density       = vis.contour;      % Draw contour lines
    plots.density.range = true;              % manually set ranges
    plots.density.x_min = 1.2;
    plots.density.x_max = 6;
    plots.density.y_min = 1.2;
    plots.density.y_max = 6;
    
    plots.density.cnt_nlev = [50 15];        % number of contour lines
    plots.density.scale_dvr = 6;             % factor that determines ranges of density
    % for contour lines to draw
    plots.density.pot_min = -0.35;           % finetuning for drawing of energy curves
    plots.density.pot_max = 0;
    
    % Plot expectation values
    plots.expect       = vis.expect;
    plots.expect.e_max = 0.3;         % manually set ranges of energy plot

end
