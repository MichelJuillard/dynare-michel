# 64 bit
./configure --with-matlab=/Volumes/MacPackages/Applications/MATLAB_R2009b.app MATLAB_VERSION=7.9

#32 bit < 7.5
./configure CFLAGS='-arch i386' CXXFLAGS='-arch i386' FFLAGS='-arch i386' LDFLAGS='-arch i386' --with-matlab=/Volumes/MacPackages/Applications/R2007a MATLAB_VERSION=7.4

#32 bit 7.5 & up
./configure CFLAGS='-arch i386' CXXFLAGS='-arch i386' FFLAGS='-arch i386' LDFLAGS='-arch i386' --with-matlab=/Volumes/MacPackages/Applications/MATLAB_R2009b_32bit/MATLAB_R2009b.app MATLAB_VERSION=7.9
