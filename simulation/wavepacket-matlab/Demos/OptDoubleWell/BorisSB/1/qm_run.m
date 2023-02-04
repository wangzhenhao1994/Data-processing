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

%% LvNE: Get ABNCD matrices
qm_setup('rho'); 
qm_init; 
qm_abncd;
qm_cleanup;

%% LvNE: Optimal control theory in full dimensionality
qm_setup('rho');
qm_init;
qm_optimal;
qm_cleanup;

