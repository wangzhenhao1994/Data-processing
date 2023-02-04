% Script to run all demos in all subdirectories. The idea is that we can
% automatically run all demos, produce all the output files, in particular
% movies, images, etc., and end up with a directory tree we can directly use 
% for our wiki webpages maintained at sourceforge.net. Use this to upload
%
% scp -r * username,wavepacket@web.sourceforge.net:htdocs/Demos/ 
%
% Should you want, for some reason, to run a subset of the demos, you can find
% similar scripts in each subdirectory.
%
% The complete run takes 9h40min on computer "bodden" without OCT
% After introducing classes: 9h55min without OCT (November 2017)
% Intel Xeon CPU E3-1241 @ 3.5GHz (March 2017)
%
% With the same hardware, a complete run takes 32h45min with OCT (V7, May 2022)

global info

%% Adiabatic: Textbook examples
cd FreeParticle;        run_demos;    cd ..;
cd LinearRamp;          run_demos;    cd ..;
cd HarmOscillator;      run_demos;    cd ..;
cd MorseOscillator;     run_demos;    cd ..;
cd DoubleWell;          run_demos;    cd ..;
cd HenonHeiles;         run_demos;    cd ..;
cd SimplePendulum;      run_demos;    cd ..;
cd GeneralPendulum;     run_demos;    cd ..;

%% Adiabatic: Molecular phsics
cd MolVibration;        run_demos;    cd ..;
cd MolTorsion;          run_demos;    cd ..;
cd MolRotation;         run_demos;    cd ..;
cd ChemReaction;        run_demos;    cd ..;
cd ElectronDynamics;    run_demos;    cd ..;

%% Non-adiabatic: Prototypical models
cd CrossingTully;       run_demos;    cd ..;
cd CrossingBerlin;      run_demos;    cd ..;
cd ConicalInter;        run_demos;    cd ..;
cd SuperExchange;       run_demos;    cd ..;

%% Non-adiabatic: Molecular physics
cd MolElectronic;       run_demos;    cd ..;
cd FemtoChem;           run_demos;    cd ..;

%% Optimal control
if strcmpi(info.system, 'Matlab') % Octave has problems with qm_matrix
    cd OptDoubleWell;       run_demos;    cd ..;
    cd OptMorseOscillator;  run_demos;    cd ..;
end
