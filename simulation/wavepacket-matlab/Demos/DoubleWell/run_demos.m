global info

%%
cd Bound_1D/1;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../..;
%%
cd Bound_1D/2;
qm_setup('wave'); qm_init(); qm_bound(); qm_cleanup();
cd ../..;

if strcmpi(info.system, 'Matlab') % Octave has problems with qm_matrix
    %%
    cd DimReduce_1D/1;
    qm_run;
    cd ../..
    %%
    cd DimReduce_1D/2;
    qm_run;
    cd ../..
    %%
    cd DimReduce_1D/3;
    qm_run;
    cd ../..
end

%%
cd Gaussian_1D/1;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Gaussian_1D/2;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
%%
cd Gaussian_1D/3;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;

