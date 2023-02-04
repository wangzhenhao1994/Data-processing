% matches pattern NL300/L0.1 etc
% stores the result in a cell structure

dirList = glob('D:\Program\simulation\H2potential\discurves\*.dat');

% Loop through the elements of the cell

for i = 1:length(dirList)

  dirname = dirList{i,1};
  [split,~,~]=schroedinger1D(3,dirname,5000,1,2.01533/2,2.01533/2);
  if split != 0
    [fPath, fName, fExt] = fileparts(dirname);
    fn =  fullfile('D:\Program\simulation\dE_hydrogen',strcat(fName,'_dE.txt'));
    dlmwrite(fn,split,'delimiter','\n')
  end
end
