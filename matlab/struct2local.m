function struct2local(S),


vnam = fieldnames(S);

for j=1:length(vnam),
  assignin('caller',vnam{j},getfield(S,vnam{j}));
end
