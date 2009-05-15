function fMessageStatus(prtfrc, njob, waitbarString, waitbarTitle, Local, MasterName, DyMo)

global funcName

if nargin<4,
  Local=1;
end

save(['comp_status_',funcName,int2str(njob)],'prtfrc','njob','waitbarString','waitbarTitle');
if Local==0,
  if isunix,
%     copyfile(['comp_status_',funcName,int2str(njob),'.mat'],['\\',MasterName,'\',DyMo(1),'$\',DyMo(4:end),'\']);
  else
    copyfile(['comp_status_',funcName,int2str(njob),'.mat'],['\\',MasterName,'\',DyMo(1),'$\',DyMo(4:end),'\']);
  end
end
