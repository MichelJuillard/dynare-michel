function y0 = getIRFRuns(jxj,jexo,type)
global M_

if nargin<3, type='metropolis'; end,

if strcmpi(type,'posterior')
  dirname = [CheckPath('metropolis') '/' ];
elseif strcmpi(type,'gsa')
  dirname = [CheckPath('GSA') '/' ];
else
  dirname = [CheckPath('prior') '/' ];
end  
filIRF = dir([dirname '/' M_.fname '_IRFs*.mat']);

y0=[];
for file = 1:length(filIRF),
  load([dirname '/' M_.fname '_IRFs' int2str(file)]);
  y0 = [y0; squeeze(STOCK_IRF(:,jxj,jexo,:))];
end
