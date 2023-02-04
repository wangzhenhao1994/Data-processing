%%
cd Bound_1D/1;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../..;

%%
cd Bound_2D/1;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../..;

%% squeezed state
cd Gaussian_1D/1;
qm_setup('wave'     ); qm_init(); qm_propa(); qm_cleanup();   
qm_setup('traj',2000); qm_init(); qm_propa(); qm_cleanup();   
cd ../..;
%% squeezed state
cd Gaussian_1D/2;
qm_setup('wave'     ); qm_init(); qm_propa(); qm_cleanup();   
qm_setup('traj',2000); qm_init(); qm_propa(); qm_cleanup();   
cd ../..;
%% coherent state
cd Gaussian_1D/3;
qm_setup('wave'     ); qm_init(); qm_propa(); qm_cleanup();   
qm_setup('traj',2000); qm_init(); qm_propa(); qm_cleanup();   
cd ../..;

%%
cd Gaussian_2D/1;
qm_setup('wave');      qm_init(); qm_propa(); qm_cleanup();   
qm_setup('traj',2000); qm_init(); qm_propa(); qm_cleanup();
cd ../..;

%%
cd Gaussian_3D/1;
qm_setup('wave');      qm_init(); qm_propa(); qm_cleanup();   
cd ../..;
