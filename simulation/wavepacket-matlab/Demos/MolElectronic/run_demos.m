global info

%% HCl+ radical photo-excitation
cd HCl+/0.5fs/BareState;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../../..;
cd HCl+/0.5fs/Floquet07;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../../..;
cd HCl+/1.0fs/BareState;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../../..;
cd HCl+/1.0fs/Floquet07;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../../..;
cd HCl+/3.0fs/BareState;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../../..;
cd HCl+/3.0fs/Floquet07;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../../..;

%% OH radical: Bare state versus Floquet dressed state representation
cd OH/0-1/BareState;   
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();   
cd ../../..;
cd OH/0-1/Dressed02;   
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();   
cd ../../..;

%% Pyrazine internal conversin: 2-state, 3-mode model
if strcmpi(info.system, 'Matlab') % Octave has problems with reduced density plots
    cd Pyrazine/1;
    qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
    cd ../..;
    cd Pyrazine/2;
    qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
    cd ../..;
end

%%
cd Retinal/1;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
cd Retinal/2;
qm_setup('wave'); qm_init(); qm_propa(); qm_cleanup();
cd ../..;
