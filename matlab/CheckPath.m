function DirectoryName = CheckPath(type)
% 06-03-2005
global M_

DirectoryName = [ M_.dname '/' type ];

if ~isdir(DirectoryName) 
    mkdir('.',DirectoryName); 
end