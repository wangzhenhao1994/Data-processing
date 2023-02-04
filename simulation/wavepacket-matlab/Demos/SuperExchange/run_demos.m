cd dia;
qm_setup('wave'); qm_init(5); qm_propa; qm_cleanup;  % Fully quantum mechanical
qm_setup('fssh'); qm_init(5); qm_propa; qm_cleanup;  % Fewest switches surface hopping
cd ..;
