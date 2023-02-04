%%
cd Bound_1D/1
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../..
%%
cd Bound_1D/2
qm_run
cd ../..

%% Short-time simulation of dephasing: Strang
cd Gaussian_1D/1;
qm_setup('wave'     ); qm_init(); qm_propa(); qm_cleanup();   
qm_setup('traj',2000); qm_init(); qm_propa(); qm_cleanup();   
cd ../..;
%% Long-time simulation of revivals: Chebychev
cd Gaussian_1D/2;
qm_setup('wave');      qm_init(); qm_propa('cheby_real',0,1e-8); qm_cleanup();   
qm_setup('traj',2000); qm_init(); qm_propa();                    qm_cleanup();   
cd ../..;
