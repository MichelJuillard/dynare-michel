function [z]=ff2a_(y,dr)
z=[];
// Copyright (C) 2001 Michel Juillard
// 
z = fff_(y)+dr('delta_s');
