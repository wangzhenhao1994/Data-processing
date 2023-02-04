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

%% LvNE: Balancing transformation
qm_setup('rho');
qm_init;
qm_balance; 
qm_cleanup;

%% LvNE: Perform truncation
qm_setup('rho');
qm_init;
qm_truncate('t',170);
qm_cleanup;

%% LvNE: Optimal control theory in reduced dimensionality
qm_setup('rho');
qm_init;
global state
state.save_suffix = 't170';
qm_optimal;
qm_cleanup;