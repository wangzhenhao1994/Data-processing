global info

%%
cd Cl-NH3/1;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Cl-NH3/2;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;

%%
cd H3+;
qm_run;
cd ..

%%
cd HF/1;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd HF/2;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;

%%
cd OH/1;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
if strcmpi(info.system, 'Matlab')
    %%
    cd OH/2;
    qm_run;
    cd ../..;
    %%
    cd OH/3;
    qm_run;
    cd ../..;
    %%
    cd OH/4;
    qm_run;
    cd ../..;
end
