## Copyright (C) 2009 Dynare Team
##
## This file is part of Dynare.
##
## Dynare is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Dynare is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

## First argument is MOD file (possibly with a path prefix)
## Second argument is path to Dynare installation to be checked
## Third argument is Dynare version to be checked

## Ask gnuplot to create graphics in text mode
## Note that setenv() was introduced in Octave 3.0.2, for compatibility
## with MATLAB
putenv("GNUTERM", "dumb")

[directory, name, ext] = fileparts(argv(){1});

printf("TEST: %s...\n", name)

addpath(argv(){2})

if !strcmp(dynare_version(), argv(){3})
  error("Incorrect version of Dynare is being tested")
endif

cd(directory)

dynare(name)

## Local variables:
## mode: Octave
## End:
