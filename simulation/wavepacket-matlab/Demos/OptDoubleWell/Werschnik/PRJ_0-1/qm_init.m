function qm_init
global hamilt space plots state time control 

prt.disp ( '***********************************' )
prt.disp ( 'Optimally driven double-well system' )
prt.disp ( '                                   ' )
prt.disp ( 'see J. Werschnik and E. K. U. Gross' )
prt.disp ( 'J. Phys. B 40, R175  (2007)        ' )
prt.disp ( 'doi:10.1088/0953-4075/40/18/R01    ' )
prt.disp ( '***********************************' )

% Save wavefunctions etc to data files
state.save_export = true; 

% Spatial discretization
space.dof{1}       = dof.fft;            % using FFT grid
space.dof{1}.mass  = 1;                  % mass
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max = +10.0;              % Upper bound of grid

% Hamiltonian operator 
hamilt.pot{1,1}        = pot.taylor;     % Quartic double well potential
hamilt.pot{1,1}.coeffs = [0; -2/4; +6/256; +24/64]; 

hamilt.dip{1}{1,1}        = dip.taylor;  % Dipole moment: Taylor series
hamilt.dip{1}{1,1}.coeffs = 1;           % Linear dipole moment, slope 1

hamilt.eigen.cutoff =   0.0;             % Cut-off entries of Hamiltonian matrix
hamilt.truncate.e_min = -05.0;           % Lower truncation of energy
hamilt.truncate.e_max = +20.0;           % Upper truncation of energy

% Calculate and save (bound) eigen states (==> qm_bound)
hamilt.eigen.start = 0;                  % Lower index
hamilt.eigen.stop  = 20;                 % Upper index

% Temporal discretization
time.steps.m_start = 000;                % Index of initial time step
time.steps.m_stop  = 040;                % Index of final time step
time.steps.m_delta = 010;                % Size of main time steps 
time.steps.s_number = 050;               % Number of substeps per main step

% Electric field pulses: Initial guess ==> qm_optimal
% Works best with ampli<0 (otherwise fluence too high)
time.pulse{1}       = efi.recta;         % Shape of envelope
time.pulse{1}.delay = 200;               % Time delay of pulse center
time.pulse{1}.fwhm  = 400;               % Pulse length
time.pulse{1}.ampli = -0.05;             % Amplitude of electric field         

% frequency resolved optical gating
time.frog.choice    = 'Gauss';           % choice of methods
time.frog.zoom      = 15;                % zoom factor for frequency axis
time.frog.width     = 20;                % zooming frequency axis

% Initial density
time.ket.choice = 'pure';                % starting from a pure state
time.ket.pure = 0;                       % selecting ground state

% Define output observables and choose control targets
control.observe.types = 'prj';           % choose types of observables
control.observe.choices={0 1 2:20}; 
control.observe.labels={'Left well' 'Right well' 'Delocalized'};
control.observe.targets=1:3; 

% Optimal control theory
control.optimal.terminal = 2;            % Which observable to be optimized
control.optimal.tolerance = 1e-4;        % Threshold terminating iteration
control.optimal.max_iter = 050;          % Max. number of iteration steps
control.optimal.alpha = 2.20;            % Penalty factor for laser fluence
control.optimal.eta  = 0.90;             % Calculate optimal backward field
control.optimal.zeta = 0.10;             % Calculate optimal forward field
control.optimal.order = 2;               % Error order for optimal fields

% Plot settings
control.plot.uxy = true;                 % Plot u(t), x(t), y(t)
control.plot.j12 = true;                 % Plot j_1(t), j_2(t), and total
control.plot.psd = true;                 % Plot power spectral density
control.plot.mov = true;                 % Animation of u(t), x(t), y(t)

% Plot u(t), x(t), y(t)
plots.control = vis.uxy;