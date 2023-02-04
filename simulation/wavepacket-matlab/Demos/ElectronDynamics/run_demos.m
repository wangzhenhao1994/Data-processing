%% 
cd PairedQuantumDots/Bound_1el;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../..;

%%
cd PairedQuantumDots/Bound_2el;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../..;

%%
cd PairedQuantumDots/Scattering_2el;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
