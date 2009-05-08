function fMessageStatus(prtfrc, njob, waitbarString, waitbarTitle, Local, MasterName, DyMo)

global fname

if nargin<4,
  Local=1;
end

save(['comp_status_',fname,int2str(njob)],'prtfrc','njob','waitbarString','waitbarTitle');
if Local==0,
  copyfile(['comp_status_',fname,int2str(njob),'.mat'],['\\',MasterName,'\',DyMo(1),'$\',DyMo(4:end),'\']);
end
