%%
cd Coherent/pi1;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Coherent/pi2;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Coherent/pi8;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;

%%
cd Double/Even;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Double/Odd;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Double/Local;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;

%%
cd Squeezed/V000;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Squeezed/V001;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Squeezed/V010;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;

%%
cd Stationary/Even;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../..;
%%
cd Stationary/Odd;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../..;
