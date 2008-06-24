function DirectoryName = CheckPath(type)

% function DirectoryName = CheckPath(type)
% Creates the repertory 'type' if it does not exist yet
%
% INPUTS
%    type
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2005-2007)
% Gnu Public License.

global M_

DirectoryName = [ M_.dname '/' type ];

if ~isdir(M_.dname)
  mkdir('.', M_.dname);
end

if ~isdir(DirectoryName)
    mkdir('.',DirectoryName);
end