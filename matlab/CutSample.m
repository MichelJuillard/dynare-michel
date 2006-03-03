function CutSample()
% stephane.adjemian@ens.fr [09-09-2005]
global M_ options_ estim_params_

npar = estim_params_.np+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nvx;

DirectoryName = CheckPath('metropolis');
file = dir([ DirectoryName '/'  M_.fname '_mh_history.mat']);
files = dir([ DirectoryName '/' M_.fname '_mh*.mat' ]);
if ~length(files)
  disp('MH:: FAILURE! there is no MH file to load here!')
  return
end
if ~length(file)
  disp('MH:: FAILURE! there is no MH-history file!')
  return
else
  load([ DirectoryName '/'  M_.fname '_mh_history'])
end
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8)
FirstDraw = floor(options_.mh_drop*TotalNumberOfMhDraws);
FirstMhFile = ceil(FirstDraw/MAX_nruns);
FirstLine = FirstDraw-(FirstMhFile-1)*MAX_nruns+1;
record.KeepedDraws.FirstMhFile = FirstMhFile;
record.KeepedDraws.FirstLine = FirstLine;
if (TotalNumberOfMhFiles-1)-(FirstMhFile+1)+1 > 0
  record.KeepedDraws.Distribution = [ MAX_nruns-FirstLine+1 ; ...
		    ones((TotalNumberOfMhFiles-1)-(FirstMhFile+1)+1,1)*MAX_nruns ; ...
		    record.MhDraws(end,3) ];
elseif TotalNumberOfMhFiles == 1
  record.KeepedDraws.Distribution = [];
elseif TotalNumberOfMhFiles == 2 & FirstMhFile > 1
  record.KeepedDraws.Distribution = [MAX_nruns-FirstLine+1 ; record.MhDraws(end,3)];  
end
save([DirectoryName '/' M_.fname '_mh_history'],'record');
fprintf('MH: Total number of Mh draws: %d.\n',TotalNumberOfMhDraws);
fprintf('MH: Total number of generated Mh files: %d.\n',TotalNumberOfMhFiles);
fprintf('MH: I''ll use mh-files %d to %d.\n',FirstMhFile,TotalNumberOfMhFiles);
fprintf('MH: In mh-file number %d i''ll start at line %d.\n',FirstMhFile,FirstLine);
fprintf('MH: Finally I keep %d draws.\n',TotalNumberOfMhDraws-FirstDraw);
disp(' ');