%%
cd C9A/1;
qm_run;
cd ../..;
%%
cd C9A/2;
qm_run;
cd ../..;
%%
cd C9A/3;
qm_run;
cd ../..;

%%
cd Fulvene/1;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../..;
%%
cd Fulvene/2;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../..;
%%
cd Fulvene/3;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Fulvene/4;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
