%------------------------------------------------------------------------------
%
% Restarts an interrupted calculation at a given time step.
%
% If you run a calculation and save the results, and if this calculation is
% interrupted unexpectedly, you can use this script to restart the calculation.
% After you have regenerated some global information using the regenerate()
% script, this script starts the calculation from some time step that was
% already saved.
%
% Input parameters are the directory where the saved calculations reside
% (wave.save.dir of the interrupted calculation), the file name of the saved files
% (wave.save.file), and the step from which you restart the calculation.
%
% Note that a handful of variables will be broken, since we leave out all steps
% in-between. These are all of expect and uncertain (they have gaps where we
% skipped the calculation of expectation values), time.steps.acf and
% time.efield.grid.* (wrong indices), time.sub.freq and time.sub.spec
% (inconsistent data)
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt
%               2007-2012 Ulf Lorenz
%
% see the README file for license details.

function qm_rerun(dir, file, start_step)

global time;

% Initializes general information and sets up log files.
prt.init (mfilename('fullpath'));

prt.disp (' ')
prt.disp ('-------------------------------------------------------------')
prt.disp (' Numerical solution of coupled time-dependent Schr�dinger    ')
prt.disp (' equations using atomic units, i.e. hbar = m_e = e = 1       ')
prt.disp (' for quantum systems interacting with electrical fields      ')
prt.disp (' (using semiclassical dipole approximation)                  ')
prt.disp (' using partial differential equations (pde) solvers          ')
prt.disp ('                                                             ')
prt.disp ('    d                 (        d    )                        ')
prt.disp (' i -- psi( R, t ) = H ( R, -i --, t ) psi( R, t )            ')
prt.disp ('   dt                 (       dR    )                        ')
prt.disp ('                                                             ')
prt.disp (' with psi(R,t)           Wavefunction      (vector, complex) ')
prt.disp (' with H = T + V -iW -F*D Hamilton operator                   ')
prt.disp (' with T = T(-i d/dR)     Kinetic energy    (scalar)          ')
prt.disp (' with V = V(R)           Potential energy  (matrix)          ')
prt.disp (' with W = W(R)           Negative imag. potential (scalar)   ')
prt.disp (' with D = D(R)           Dipole moment     (matrix)          ')
prt.disp (' with F = F(t)           Electrical field  (scalar, time-dep)')
prt.disp ('                                                             ')
prt.disp ('-------------------------------------------------------------')
prt.disp (' ')


%% Initialize coupling scheme and spatial/temporal discretization
init.grid;

% Initialize the electric field
init.efield;

% Initialize Hamiltonian operator 
init.hamilt;

% Initialize temporal discretization
tmp.init;

% Initialize wave function
wave.init;

% Initialize expectation values and uncertainties of observables
obs.init;

%% Run one step for the setup
for step = 1
    
    % Numerical propagation using pde solvers, possibly with absorbing boundary conditions
    wave.propa ( step );

    % Transform to adiabatic representation (if desired)
    wave.adiabatic ( 'dia2adi' );
    
    % Expectation values and uncertainties
    wave.observe ( step );
    
    % Logging and plot title
    obs.log ( step );       

    % Get spectrum as Fourier transform of autocorrelation
    wave.spectrum ( step );

    % Plot densities and expectation values
    plot.ket ( step );

    % Transform back to diabatic representation
    wave.adiabatic ( 'adi2dia' );

    % Store the wave function.
    wave.save ( step );
    
% End of main loop
end


% Load the wave function
context = wave.load(dir, file);
wave.load(context, start_step, true);
clear context


% And continue from where we left off.
for step = start_step+1 : time.steps.m_number
    
    % Numerical propagation using pde solvers, possibly with absorbing boundary conditions
    wave.propa ( step );

    % Transform to adiabatic representation (if desired)
    wave.adiabatic ( 'dia2adi' );
    
    % Expectation values and uncertainties
    wave.observe ( step );
    
    % Logging and plot title
    obs.log ( step );       

    % Get spectrum as Fourier transform of autocorrelation
    wave.spectrum ( step );

    % Plot densities and expectation values
    plot.ket ( step );

    % Transform back to diabatic representation
    wave.adiabatic ( 'adi2dia' );

    % Store the wave function.
    wave.save ( step );
    
% End of main loop
end

% Output clock/date/time
prt.clock;

end