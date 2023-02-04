%% Solve TISE
qm_setup('wave');
qm_init; 
qm_bound;            
qm_cleanup;

%% Matrix representations
qm_setup('wave'); 
qm_init; 
qm_matrix; 
qm_cleanup;

%% TDSE: Get ABNCD matrices
qm_setup('rho'); 
qm_init; 
qm_abncd;
qm_cleanup;

%% TDSE: Optimal Control Theory
qm_setup('rho');
qm_init;
qm_optimal;
qm_cleanup;
