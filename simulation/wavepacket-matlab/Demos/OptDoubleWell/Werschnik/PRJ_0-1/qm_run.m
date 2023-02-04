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

%% Switch to ABNCD control language
qm_setup('ket');  
qm_init; 
qm_abncd;      
qm_cleanup;

%% Optimal control
qm_setup('ket');  
qm_init; 
qm_optimal; 
qm_cleanup;