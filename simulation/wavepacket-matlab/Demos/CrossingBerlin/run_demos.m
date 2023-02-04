%% Example published in 2002
cd Double/1;
qm_setup('fssh',  3000); qm_init(); qm_propa(); qm_cleanup(); % Tully's fewest switches  
qm_setup('sssh',2,3000); qm_init(); qm_propa(); qm_cleanup(); % Landau-Zener variant 2 
qm_setup('wave'       ); qm_init(); qm_propa(); qm_cleanup(); % Quantum wave packet  
cd ../..;
%% Example published in 2004
cd Double/2;
qm_setup('fssh',  3000); qm_init(); qm_propa(); qm_cleanup();   
qm_setup('sssh',2,3000); qm_init(); qm_propa(); qm_cleanup();   
qm_setup('wave'       ); qm_init(); qm_propa(); qm_cleanup();   
cd ../..;


%% Propagations
cd Single/1;
qm_setup('fssh',  3000); qm_init(); qm_propa(); qm_cleanup(); % Tully's fewest switches
qm_setup('sssh',2,3000); qm_init(); qm_propa(); qm_cleanup(); % Landau-Zener variant 2
qm_setup('wave'       ); qm_init(); qm_propa(); qm_cleanup(); % Quantum wave packet
cd ../..;

%% Bound state calculations
cd Single/2;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();   
cd ../..;
