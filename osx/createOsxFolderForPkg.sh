#!/bin/bash

set -ex

TOP_DIR=/Users/Houtan/Documents/DYNARE
TOP_DYN_DIR=$TOP_DIR/dynare

VERSION=4.2.2
INSTALLDIRNAME=dynare-$VERSION-osx
INSTALLDIR=$TOP_DIR/$INSTALLDIRNAME
mkdir $INSTALLDIR

PATH=${PATH}:/Applications/Octave.app/Contents/Resources/bin:/usr/local/bin:/usr/local/sbin:/usr/texbin:/usr/local/include:/usr/local/lib:/usr/local/libexec:/usr/local/share:/usr/local/Cellar/

########################
# UPDATE DYNARE SOURCE #
########################
cd $TOP_DYN_DIR
autoreconf -si


########################
# BEGIN MAKING PACKAGE #
########################
# create directories
mkdir "$INSTALLDIR/doc"
mkdir "$INSTALLDIR/doc/dynare++"
mkdir "$INSTALLDIR/doc/dynare.html"
mkdir "$INSTALLDIR/dynare++"
mkdir "$INSTALLDIR/mex"
mkdir "$INSTALLDIR/mex/octave"
mkdir "$INSTALLDIR/mex/matlab"
mkdir "$INSTALLDIR/mex/matlab/osx64"
mkdir "$INSTALLDIR/mex/matlab/osx32-7.4"
mkdir "$INSTALLDIR/mex/matlab/osx32-7.5-7.13"

# top level
cp $TOP_DYN_DIR/dynare.el                                        $INSTALLDIR
cp $TOP_DYN_DIR/license.txt                                      $INSTALLDIR

# matlab
cp -r $TOP_DYN_DIR/matlab                                        $INSTALLDIR

# examples
cp -r $TOP_DYN_DIR/examples                                      $INSTALLDIR


##########################################################
# FIRST BUILD 32 BIT EVERYTHING, 32 BIT MATLAB < 7.5 MEX #
##########################################################
./configure CFLAGS='-arch i386' CXXFLAGS='-arch i386' FFLAGS='-arch i386' LDFLAGS='-arch i386' --with-matlab=/Applications/MATLAB/R2007a MATLAB_VERSION=7.4 --with-gsl=/usr/local/Cellar/gsl/1.15_32bit
make pdf
make


########################
# MAKE BULK OF PACKAGE #
########################
# compiled preprocessor
cp $TOP_DYN_DIR/preprocessor/dynare_m                            $INSTALLDIR/matlab

# Matlab
cp $TOP_DYN_DIR/mex/build/matlab/block_kalman_filter/*.mexmaci   $INSTALLDIR/mex/matlab/osx32-7.4
cp $TOP_DYN_DIR/mex/build/matlab/bytecode/*.mexmaci              $INSTALLDIR/mex/matlab/osx32-7.4
cp $TOP_DYN_DIR/mex/build/matlab/dynare_simul_/*.mexmaci         $INSTALLDIR/mex/matlab/osx32-7.4
cp $TOP_DYN_DIR/mex/build/matlab/estimation/*.mexmaci            $INSTALLDIR/mex/matlab/osx32-7.4
cp $TOP_DYN_DIR/mex/build/matlab/gensylv/*.mexmaci               $INSTALLDIR/mex/matlab/osx32-7.4
cp $TOP_DYN_DIR/mex/build/matlab/k_order_perturbation/*.mexmaci  $INSTALLDIR/mex/matlab/osx32-7.4
cp $TOP_DYN_DIR/mex/build/matlab/kalman_steady_state/*.mexmaci   $INSTALLDIR/mex/matlab/osx32-7.4
cp $TOP_DYN_DIR/mex/build/matlab/kronecker/*.mexmaci             $INSTALLDIR/mex/matlab/osx32-7.4
cp $TOP_DYN_DIR/mex/build/matlab/mjdgges/*.mexmaci               $INSTALLDIR/mex/matlab/osx32-7.4
cp $TOP_DYN_DIR/mex/build/matlab/ms_sbvar/*.mexmaci              $INSTALLDIR/mex/matlab/osx32-7.4

# Octave
cp $TOP_DYN_DIR/mex/build/octave/block_kalman_filter/*.mex       $INSTALLDIR/mex/octave
cp $TOP_DYN_DIR/mex/build/octave/bytecode/*.mex                  $INSTALLDIR/mex/octave
cp $TOP_DYN_DIR/mex/build/octave/dynare_simul_/*.mex             $INSTALLDIR/mex/octave
cp $TOP_DYN_DIR/mex/build/octave/estimation/*.mex                $INSTALLDIR/mex/octave
cp $TOP_DYN_DIR/mex/build/octave/gensylv/*.mex                   $INSTALLDIR/mex/octave
cp $TOP_DYN_DIR/mex/build/octave/k_order_perturbation/*.mex      $INSTALLDIR/mex/octave
cp $TOP_DYN_DIR/mex/build/octave/kalman_steady_state/*.mex       $INSTALLDIR/mex/octave
cp $TOP_DYN_DIR/mex/build/octave/kronecker/*.mex                 $INSTALLDIR/mex/octave
cp $TOP_DYN_DIR/mex/build/octave/mjdgges/*.mex                   $INSTALLDIR/mex/octave
cp $TOP_DYN_DIR/mex/build/octave/ms_sbvar/*.mex                  $INSTALLDIR/mex/octave
cp $TOP_DYN_DIR/mex/build/octave/ordschur/*.oct                  $INSTALLDIR/mex/octave
cp $TOP_DYN_DIR/mex/build/octave/qzcomplex/*.oct                 $INSTALLDIR/mex/octave

# dynare++
cp $TOP_DYN_DIR/dynare++/src/dynare++                            $INSTALLDIR/dynare++
cp $TOP_DYN_DIR/dynare++/extern/matlab/dynare_simul.m            $INSTALLDIR/dynare++

# doc
cp $TOP_DYN_DIR/doc/bvar-a-la-sims.pdf                           $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/dr.pdf                                       $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/dynare.pdf                                   $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/guide.pdf                                    $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/macroprocessor/macroprocessor.pdf            $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/parallel/parallel.pdf                        $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/preprocessor/preprocessor.pdf                $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/userguide/UserGuide.pdf                      $INSTALLDIR/doc

# doc (dynare++)
cp $TOP_DYN_DIR/dynare++/doc/dynare++-tutorial.pdf               $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/doc/dynare++-ramsey.pdf                 $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/sylv/sylvester.pdf                      $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/tl/cc/tl.pdf                            $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/integ/cc/integ.pdf                      $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/kord/kord.pdf                           $INSTALLDIR/doc/dynare++


##############################################
# RETURN TO BUILD 32 BIT MATLAB 7.5 & UP MEX #
##############################################
make clean
cd $TOP_DYN_DIR/mex/build/matlab
./configure CFLAGS='-arch i386' CXXFLAGS='-arch i386' FFLAGS='-arch i386' LDFLAGS='-arch i386' --with-matlab=/Applications/MATLAB/MATLAB_R2009b_32bit/MATLAB_R2009b.app MATLAB_VERSION=7.9 MEXEXT='mexmaci' --with-gsl=/usr/local/Cellar/gsl/1.15_32bit
make

# Matlab
cp $TOP_DYN_DIR/mex/build/matlab/block_kalman_filter/*.mexmaci   $INSTALLDIR/mex/matlab/osx32-7.5-7.13
cp $TOP_DYN_DIR/mex/build/matlab/bytecode/*.mexmaci              $INSTALLDIR/mex/matlab/osx32-7.5-7.13
cp $TOP_DYN_DIR/mex/build/matlab/dynare_simul_/*.mexmaci         $INSTALLDIR/mex/matlab/osx32-7.5-7.13
cp $TOP_DYN_DIR/mex/build/matlab/estimation/*.mexmaci            $INSTALLDIR/mex/matlab/osx32-7.5-7.13
cp $TOP_DYN_DIR/mex/build/matlab/gensylv/*.mexmaci               $INSTALLDIR/mex/matlab/osx32-7.5-7.13
cp $TOP_DYN_DIR/mex/build/matlab/k_order_perturbation/*.mexmaci  $INSTALLDIR/mex/matlab/osx32-7.5-7.13
cp $TOP_DYN_DIR/mex/build/matlab/kalman_steady_state/*.mexmaci   $INSTALLDIR/mex/matlab/osx32-7.5-7.13
cp $TOP_DYN_DIR/mex/build/matlab/kronecker/*.mexmaci             $INSTALLDIR/mex/matlab/osx32-7.5-7.13
cp $TOP_DYN_DIR/mex/build/matlab/mjdgges/*.mexmaci               $INSTALLDIR/mex/matlab/osx32-7.5-7.13
cp $TOP_DYN_DIR/mex/build/matlab/ms_sbvar/*.mexmaci              $INSTALLDIR/mex/matlab/osx32-7.5-7.13


#####################################
# RETURN TO BUILD 64 BIT MATLAB MEX #
#####################################
make clean
cd $TOP_DYN_DIR/mex/build/matlab
./configure --with-matlab=/Applications/MATLAB/MATLAB_R2009b.app MATLAB_VERSION=7.9 --with-gsl=/usr/local/Cellar/gsl/1.15
make

# Matlab
cp $TOP_DYN_DIR/mex/build/matlab/block_kalman_filter/*.mexmaci64   $INSTALLDIR/mex/matlab/osx64
cp $TOP_DYN_DIR/mex/build/matlab/bytecode/*.mexmaci64              $INSTALLDIR/mex/matlab/osx64
cp $TOP_DYN_DIR/mex/build/matlab/dynare_simul_/*.mexmaci64         $INSTALLDIR/mex/matlab/osx64
cp $TOP_DYN_DIR/mex/build/matlab/estimation/*.mexmaci64            $INSTALLDIR/mex/matlab/osx64
cp $TOP_DYN_DIR/mex/build/matlab/gensylv/*.mexmaci64               $INSTALLDIR/mex/matlab/osx64
cp $TOP_DYN_DIR/mex/build/matlab/k_order_perturbation/*.mexmaci64  $INSTALLDIR/mex/matlab/osx64
cp $TOP_DYN_DIR/mex/build/matlab/kalman_steady_state/*.mexmaci64   $INSTALLDIR/mex/matlab/osx64
cp $TOP_DYN_DIR/mex/build/matlab/kronecker/*.mexmaci64             $INSTALLDIR/mex/matlab/osx64
cp $TOP_DYN_DIR/mex/build/matlab/mjdgges/*.mexmaci64               $INSTALLDIR/mex/matlab/osx64
cp $TOP_DYN_DIR/mex/build/matlab/ms_sbvar/*.mexmaci64              $INSTALLDIR/mex/matlab/osx64

# clean everything
cd $TOP_DYN_DIR
make distclean

# remove .DS_Store files
cd $INSTALLDIR
find . -name *.DS_Store -type f -exec rm {} \;

# adjust permissions
cd $TOP_DIR
chmod -R g+w $INSTALLDIR

echo "DONE :)"
