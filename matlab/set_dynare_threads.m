function not = set_dynare_threads(n)
% This function sets the number of threads used by dynare's mex files.
%
% INPUTS 
%  o n    [integer]   scalar specifying the number of threads to be used.    
%
% OUTPUTS 
%  o not  [integer]   scalar, number of threads.    
%
% REMARKS The default value of n is the number of processors on the platform.    
    
% Copyright (C) 2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

    not = 1;
    if ~isopenmp()% This version of Dynare does not use multithreaded mex files!
        return
    end
    MaxNumberOfThreads = maxNumCompThreads();
    if ~nargin% Default.
        not = MaxNumberOfThreads;
        setenv('DYNARE_NUM_THREADS',int2str(MaxNumberOfThreads));
    else
        if (n>MaxNumberOfThreads)
            disp(['You want to use ' int2str(n) ' threads but your platform has only ' int2str(MaxNumberOfThreads) ' processors!'])
            reply = input(['Do you really want to use ' int2str(n) ' threads ?  Yes/[No]: '],'s');
            if isempty(reply)
                reply = 'No';
            end
            if strcmpi(reply,'No')
                nn = input(['Choose a number of threads between 1 and [' int2str(MaxNumberOfThreads) ']: ']);
                if isempty(nn)
                    nn = MaxNumberOfThreads;
                end
                if (nn>MaxNumberOfThreads)
                    disp(['To my knowledge ' int2str(nn) ' is greater than ' int2str(MaxNumberOfThreads) '!...'])
                    disp(' ')
                    not = set_dynare_threads(n);
                    return
                end
                not = nn;
                setenv('DYNARE_NUM_THREADS',int2str(nn));
            elseif strcmpi(reply,'Yes')
                not = n;
                setenv('DYNARE_NUM_THREADS',int2str(n));
            else
                disp(['You have to answer by Yes or No...'])
                disp(' ')
                not = set_dynare_threads(n);
                return
            end
        else
            not = n;
            setenv('DYNARE_NUM_THREADS',int2str(n)); 
        end
    end