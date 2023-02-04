%% Solve TISE in position/momentum representation 
qm_setup('wave'); 
qm_init(); 
qm_bound(); 
qm_cleanup();  

%% Get matrix elements
qm_setup('wave');  
qm_init(); 
qm_matrix(); 
qm_cleanup();

%% Transform to energy representation
qm_setup('rho'); 
qm_init(); 
qm_abncd(); 
qm_cleanup(); 

%% Solve TDSE in energy representation 
qm_setup('rho'); 
qm_init(); 
qm_propa(); 
qm_cleanup();  