#!/bin/bash

#
# BEGIN EDIT
#
DYNAREV=4.2.2
TOPDIR=../dynare-4.2.2
#
# END EDIT
#

INSTALDIR=dynare-$DYNAREV
mkdir $INSTALDIR


#
# TOP LEVEL
#
cp $TOPDIR/dynare.el    $INSTALDIR
cp $TOPDIR/license.txt  $INSTALDIR


#
# MATLAB
#
cp -r $TOPDIR/matlab                      $INSTALDIR
cp    $TOPDIR/preprocessor/dynare_m       $INSTALDIR/matlab


#
# MEX
#
mkdir "$INSTALDIR/mex"

# Matlab
cp -r $TOPDIR/mex/matlab              $INSTALDIR/mex

# Octave
mkdir "$INSTALDIR/mex/octave"
cp $TOPDIR/mex/build/octave/bytecode/*.mex              $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/dynare_simul_/*.mex         $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/gensylv/*.mex               $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/k_order_perturbation/*.mex  $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/kronecker/*.mex             $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/mjdgges/*.mex               $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/ordschur/*.oct              $INSTALDIR/mex/octave
cp $TOPDIR/mex/build/octave/qzcomplex/*.oct             $INSTALDIR/mex/octave


#
# EXAMPLES
#
cp -r $TOPDIR/examples $INSTALDIR/


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
cp $TOPDIR/doc/bvar-a-la-sims.pdf                 $INSTALDIR/doc
cp $TOPDIR/doc/dr.pdf                             $INSTALDIR/doc
cp $TOPDIR/doc/dynare.pdf                         $INSTALDIR/doc
cp $TOPDIR/doc/guide.pdf                          $INSTALDIR/doc
cp $TOPDIR/doc/macroprocessor/macroprocessor.pdf  $INSTALDIR/doc
cp $TOPDIR/doc/parallel/parallel.pdf              $INSTALDIR/doc
cp $TOPDIR/doc/preprocessor/preprocessor.pdf      $INSTALDIR/doc
cp $TOPDIR/doc/userguide/UserGuide.pdf            $INSTALDIR/doc

# html
mkdir "$INSTALDIR/doc/dynare.html"
cp -r $TOPDIR/doc/dynare.html/*.png               $INSTALDIR/doc/dynare.html
cp -r $TOPDIR/doc/dynare.html/*.html              $INSTALDIR/doc/dynare.html

# pdf (dynare++)
mkdir "$INSTALDIR/doc/dynare++"
cp $TOPDIR/dynare++/doc/dynare++-tutorial.pdf     $INSTALDIR/doc/dynare++
cp $TOPDIR/dynare++/doc/dynare++-ramsey.pdf       $INSTALDIR/doc/dynare++
cp $TOPDIR/dynare++/sylv/sylvester.pdf            $INSTALDIR/doc/dynare++
cp $TOPDIR/dynare++/tl/cc/tl.pdf                  $INSTALDIR/doc/dynare++
cp $TOPDIR/dynare++/integ/cc/integ.pdf            $INSTALDIR/doc/dynare++
cp $TOPDIR/dynare++/kord/kord.pdf                 $INSTALDIR/doc/dynare++

./removeDsStore.sh