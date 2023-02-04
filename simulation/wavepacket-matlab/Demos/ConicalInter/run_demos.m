%%
cd GeoPhase/1;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd GeoPhase/2;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;

%%
cd Linear/1;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Linear/2;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;

%%
cd OneDim/1;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;

%%
cd SpinBoson/1;
qm_setup('wave'); qm_init(); qm_propa('cheby_real'); qm_cleanup();
cd ../..;
%%
cd SpinBoson/2;
qm_setup('wave'); qm_init(); qm_propa('cheby_real'); qm_cleanup();
cd ../..;

%%
cd Warping/1;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Warping/2;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
