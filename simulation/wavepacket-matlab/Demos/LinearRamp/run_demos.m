cd Gaussian_1D/1;   
qm_setup('wave');      qm_init(); qm_propa(); qm_cleanup();   
qm_setup('traj',2000); qm_init(); qm_propa(); qm_cleanup();   
cd ../..;
