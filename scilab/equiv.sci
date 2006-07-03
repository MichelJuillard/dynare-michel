function [y]=equiv(x)
y=[];
// Copyright (C) 2001 Michel Juillard
// 
ys_ = x;
 
//!! Unknown function taylor ,the original calling sequence is used
taylor();
 
//!! Unknown function newt ,the original calling sequence is used
y = newt('ff20_',ys_);
