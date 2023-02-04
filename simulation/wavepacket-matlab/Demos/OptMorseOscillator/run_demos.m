%%
cd AMO_Gauss/1;
qm_run;
cd ../..; 
%%
copyfile ('AMO_Gauss/1/ket_optimal_40.dat','AMO_Gauss/2/')
cd AMO_Gauss/2;
qm_run;
cd ../..; 

%%
cd OVL_0-1/1;
qm_run;
cd ../..; 
%%
cd OVL_0-1/2;
qm_run;
cd ../..; 

%%
cd PRJ_0-1/1;
qm_run;
cd ../..; 
%%
cd PRJ_0-1/2;
qm_run;
cd ../..; 
%%
cd PRJ_0-1/3;
qm_run;
cd ../..; 

