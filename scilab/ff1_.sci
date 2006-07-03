function [y]=ff1_(x)
// Copyright (C) 2001 Michel Juillard
// 
global ex_

n1 = length(x)-exo_nbr;
ex_(it_+xkmin_-ykmin_,:) = x(n1+1:$)';
y = evstr(fname_+'_ff(x(1:n1))');
 
 
 
