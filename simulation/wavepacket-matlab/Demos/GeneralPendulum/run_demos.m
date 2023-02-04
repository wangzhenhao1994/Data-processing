%%
cd Planar/Stationary/1;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../../..;
%%
cd Planar/Stationary/2;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../../..;

%%
cd Spherical/Stationary/1;
qm_setup('wave'); qm_init (10,1); qm_bound(); qm_cleanup()
cd ../../..;
%%
cd Spherical/Stationary/2;
qm_setup('wave'); qm_init (10,1); qm_bound(); qm_cleanup()
cd ../../..;
