#!/bin/bash

#
# BEGIN EDIT
#
DYNAREV=4.3-unstable
TOPDIR=..
VIRGINTOPDIR=../../untouchedDynare430/dynare-4.3-unstable
DOCDIR=doc
DYNPPDIR=dynare++
#
# END EDIT
#

INSTALDIR=dynare-$DYNAREV
mkdir $INSTALDIR


#
# TOP LEVEL
#
cp build_dynare.m       $INSTALDIR
cp $TOPDIR/dynare.el    $INSTALDIR
cp $TOPDIR/license.txt  $INSTALDIR


#
# SRC
#
mkdir "$INSTALDIR/src"
mkdir "$INSTALDIR/src/preprocessor"
cp -r  boost_1_45_0                                   $INSTALDIR/src
cp -r  $VIRGINTOPDIR/dynare++                         $INSTALDIR/src
cp -r  $VIRGINTOPDIR/m4                               $INSTALDIR/src
cp -r  $VIRGINTOPDIR/mex                              $INSTALDIR/src
cp     $TOPDIR/mex/build/matlab/configure.ac          $INSTALDIR/src/mex/build/matlab
cp     $VIRGINTOPDIR/preprocessor/CodeInterpreter.hh  $INSTALDIR/src/preprocessor
rm -rf $INSTALDIR/src/mex/octave
rm -rf $INSTALDIR/src/mex/build/octave


#
# MATLAB
#
cp -r $VIRGINTOPDIR/matlab                $INSTALDIR
cp    $TOPDIR/preprocessor/dynare_m       $INSTALDIR/matlab


#
# MEX
#
mkdir "$INSTALDIR/mex"
mkdir "$INSTALDIR/mex/matlab"

# Octave
mkdir "$INSTALDIR/mex/octave"
cp $TOPDIR/mex/build/octave/bytecode/*.mex              $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/dynare_simul_/*.mex         $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/gensylv/*.mex               $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/k_order_perturbation/*.mex  $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/kalman_steady_state/*.mex   $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/kronecker/*.mex             $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/mjdgges/*.mex               $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/ordschur/*.oct              $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/qzcomplex/*.oct             $INSTALDIR/mex/octave


#
# DYNARE++
#
mkdir "$INSTALDIR/dynare++"
cp $TOPDIR/dynare++/src/dynare++                  $INSTALDIR/dynare++
cp $TOPDIR/dynare++/extern/matlab/dynare_simul.m  $INSTALDIR/dynare++


#
# DOC
#

# pdf (dynare)
mkdir "$INSTALDIR/doc"
cp $DOCDIR/bvar-a-la-sims.pdf                 $INSTALDIR/doc
cp $DOCDIR/dr.pdf                             $INSTALDIR/doc
cp $DOCDIR/dynare.pdf                         $INSTALDIR/doc
cp $DOCDIR/guide.pdf                          $INSTALDIR/doc
cp $DOCDIR/macroprocessor/macroprocessor.pdf  $INSTALDIR/doc
cp $DOCDIR/parallel/parallel.pdf              $INSTALDIR/doc
cp $DOCDIR/preprocessor/preprocessor.pdf      $INSTALDIR/doc
cp $DOCDIR/userguide/UserGuide.pdf            $INSTALDIR/doc

# html
mkdir "$INSTALDIR/doc/dynare.html"
cp -r $DOCDIR/dynare.html/*.png                      $INSTALDIR/doc/dynare.html
cp -r $DOCDIR/dynare.html/*.html                     $INSTALDIR/doc/dynare.html

# pdf (dynare++)
mkdir "$INSTALDIR/doc/dynare++"
cp $DYNPPDIR/doc/dynare++-tutorial.pdf     $INSTALDIR/doc/dynare++
cp $DYNPPDIR/doc/dynare++-ramsey.pdf       $INSTALDIR/doc/dynare++
cp $DYNPPDIR/sylv/sylvester.pdf            $INSTALDIR/doc/dynare++
cp $DYNPPDIR/tl/cc/tl.pdf                  $INSTALDIR/doc/dynare++
cp $DYNPPDIR/integ/cc/integ.pdf            $INSTALDIR/doc/dynare++
cp $DYNPPDIR/kord/kord.pdf                 $INSTALDIR/doc/dynare++

#
# EXAMPLES
#
cp -r $VIRGINTOPDIR/examples $INSTALDIR/