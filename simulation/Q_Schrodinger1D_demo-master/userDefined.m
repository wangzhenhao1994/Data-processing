%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pot = load('D:\SF_FFT\Program\simulation\H2potential\discurves\1_0_1.dat');
Vb=(pot(1,2)-min(pot(:,2)))*27.211386245988;                 % potential barrier [eV]
a=1e15;
Ltot=10*5.29177210903*1e-11;

z=0:dz:Ltot;
%V0=a*(z-Ltot/2).^2;

V0 = csapi(pot(:,1)*5.29177210903*1e-11,(pot(:,2)-min(pot(:,2)))*27.211386245988,z);
%V0(V0>Vb)=Vb;

if TM_Method==1
  display('ERROR: TM_Method not valid for that potential')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%