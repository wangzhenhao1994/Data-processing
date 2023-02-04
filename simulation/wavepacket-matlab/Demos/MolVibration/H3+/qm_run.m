% Calculates the lowest three vibrational states of H3+. See work by
% W. Meyer, P. Botschwina, P. Burton, J. Chem. Phys. 84:891 (1986)
% According to Table VIII of that reference (doi:10.1063/1.450534),
% we should see the eigenstates at 4345, 6861, 7530 cm^-1.
% Due to convergence problems, our results may deviate by about
% 10-20 cm^-1.
%
% This script uses propagation in imaginary time to get the
% eigenstates. It will first calculate the groundstate, then the 
% first excited state, then the second one. Each of the latter two
% calculations uses a the 'cheby_imag' imaginary time propagator
% that propagates only in the space perpendicular to the lower 
% (i.e. previously calculated) eigenstates. 

global atomic hamilt

% Do setup only once, in order not to loose state.bound{.}
qm_setup('wave');

% Do initialization before main loop => hamilt.eigen
qm_init(); 

% Allocate
energies = zeros(hamilt.eigen.stop-hamilt.eigen.start+1,1);

% Main loop over eigenstates
for i_eig = hamilt.eigen.start:hamilt.eigen.stop
    
    % Get bound state by propagation in imaginary time
    qm_propa('cheby_imag',0,1e-8); 
    qm_cleanup();
    
    % Retrieve and convert energies to cm^-1
    global expect
    energies(i_eig+1) = expect.total(end)  * atomic.w.cm_1;
    
end

% Display all bound state energies
prt.disp (' ')
prt.disp( 'Bound state energies (cm^-1):' )
prt.disp (' ')
for i_eig = hamilt.eigen.start:hamilt.eigen.stop
    disp( [int2str(i_eig) ': ' num2str(energies(i_eig+1))] )
end

