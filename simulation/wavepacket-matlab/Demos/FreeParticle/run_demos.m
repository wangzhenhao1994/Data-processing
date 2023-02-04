%% zero initial momentum
cd Gaussian_1D/1;   
qm_setup('wave');      qm_init(); qm_propa(); qm_cleanup();   
qm_setup('traj',2000); qm_init(); qm_propa(); qm_cleanup();   
cd ../..;
%% non-zero initial momentum
cd Gaussian_1D/2;   
qm_setup('wave');      qm_init(); qm_propa(); qm_cleanup();   
qm_setup('traj',2000); qm_init(); qm_propa(); qm_cleanup();   
cd ../..;


%%
cd Gaussian_2D/1;   
qm_setup('wave');      qm_init(); qm_propa(); qm_cleanup();
qm_setup('traj',2000); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Gaussian_2D/2;
qm_setup('wave');      qm_init(); qm_propa(); qm_cleanup();   
qm_setup('traj',2000); qm_init(); qm_propa(); qm_cleanup();
cd ../..;


%%
cd SchroediCat/1;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd SchroediCat/2;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
