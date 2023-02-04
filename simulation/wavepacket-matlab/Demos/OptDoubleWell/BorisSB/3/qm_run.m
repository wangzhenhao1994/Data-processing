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
qm_cleanup();

%% LvNE: H2 optimal model reduction
qm_setup('rho');
qm_init;
qm_H2model(170); 
qm_cleanup;

%% LvNE: Optimal control theory in reduced dimensionality
qm_setup('rho');
qm_init;
global state
state.save_suffix = 'h170';
qm_optimal;
qm_cleanup;