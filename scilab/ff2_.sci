function [y]=ff2_(x)
y=[];
// Copyright (C) 2001 Michel Juillard
// 
global('np_','endo_nbr','nf_','ex_','it_','ykmin_','xkmin_','exo_nbr','jp_','jf_','kf_','exe_');
y = zeros(np_+endo_nbr+nf_,1);
%v = x(1:np_)
y(1:np_,1) = %v(:);
ex_(it_-ykmin_+xkmin_,:) = x(np_+1:np_+exo_nbr)';
 
//!! Unknown function h1_ ,the original calling sequence is used
%v = h1_(y(1:np_),ex_(it_-ykmin_+xkmin_,:)')
y(np_+jp_,1) = %v(:);
 
//!! Unknown function g1_ ,the original calling sequence is used
%v = g1_(y(1:np_),ex_(it_-ykmin_+xkmin_,:)')
y(np_+jf_,1) = %v(:);
 
//!! Unknown function g1_ ,the original calling sequence is used
%v = g1_(y(np_+jp_),exe_)
y(kf_,1) = %v(:);
y = ff_(y);
 
 
 
