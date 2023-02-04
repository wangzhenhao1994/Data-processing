global hamilt

% Do setup only once, in order not to loose state.bound{.}
qm_setup('wave');

% Do initialization before main loop => hamilt.eigen
qm_init(false); 

% Allocate
energies = zeros(hamilt.eigen.stop-hamilt.eigen.start+1,1);

% Main loop over eigenstates
for i_eig = hamilt.eigen.start:hamilt.eigen.stop
    prt.disp ('*****************************')
    prt.disp ( ['Calculating state ', i_eig] )
    prt.disp ('*****************************')
    prt.disp (' ')
    
    % Get bound state by propagation in imaginary time
    if i_eig==hamilt.eigen.start
        qm_init(true); % with graphics
    else
        qm_init(false); % without graphics
    end
    qm_propa('cheby_imag',0,1e-8); 
    qm_cleanup();
    
    % Retrieve and convert energies to cm^-1
    global expect
    energies(i_eig+1) = expect.total(end);
    
end

% Display all bound state energies
prt.disp (' ')
prt.disp ('Bound state energies:' )
prt.disp (' ')
for i_eig = hamilt.eigen.start:hamilt.eigen.stop
    disp( [int2str(i_eig) ': ' num2str(energies(i_eig+1))] )
end


