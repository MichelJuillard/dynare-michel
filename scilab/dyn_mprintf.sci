function []=dyn_mprintf(fmt,varargin)
disp(fmt)
disp(varargin)
mprintf(fmt,varargin);
mfprintf(fh_log,fmt,varargin);
