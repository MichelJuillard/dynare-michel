function fid = openfile(filename,permission)
%Open a file in read/write mode, catching errors
%  FID = OPENFILE(FILENAME,PERMISSION) opens file FILENAME
%  in PERMISSION mode ('r' or 'w') and return a file identifier FID.

%  Copyright (C) 2003 Guillaume Flandin <Guillaume@artefact.tk>
%  $Revision: 1.1.2.1 $Date: 2004/03/31 14:46:22 $

[fid, errmsg] = fopen(filename,permission);
if ~isempty(errmsg)
	switch permission
		case 'r'
			error(sprintf('Cannot open %s in read mode.',filename));
		case 'w'
			error(sprintf('Cannot open %s in write mode.',filename));
		otherwise
			error(errmsg);
	end
end
