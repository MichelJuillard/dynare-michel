#!/bin/bash

declare -i total=0;
declare -i failed=0;
declare -a failed_tests=("");

# Parse TRS Files
for file in $1 ; do
  # Find number of tests run in trs file
  ((total += `grep number-tests $file | cut -d: -f3`))

  # Find number of tests failed in trs file
  numfailed=`grep number-failed-tests $file | cut -d: -f3`
  if [ $numfailed -ne 0 ] ; then
    ((failed += $numfailed))
    for failedfile in `grep list-of-failed-tests $file | cut -d: -f3` ; do
      failed_tests=("${failed_tests[@]}" "$failedfile");
    done
  fi
done
((passed=$total-$failed));

# Determine if we are parsing Matlab or Octave trs files
if [ `grep -c '.m.trs' <<< $1` -eq 0 ]; then
  prg='OCTAVE';
  outfile='run_test_octave_output.txt'
else
  prg='MATLAB';
  outfile='run_test_matlab_output.txt'
fi

# Print Output
echo '================================'             > $outfile
echo 'DYNARE MAKE CHECK '$prg' RESULTS'            >> $outfile
echo '================================'            >> $outfile
echo '| TOTAL: '$total                             >> $outfile
echo '|  PASS: '$passed                            >> $outfile
echo '|  FAIL: '$failed                            >> $outfile
if [ $failed -gt 0 ] ; then
  echo '|  LIST OF FAILED TESTS:'                  >> $outfile
  for file in ${failed_tests[@]} ; do
    if [ "$prg" == "MATLAB" ]; then
      modfile=`sed 's/\.m\.trs/\.mod/g' <<< $file` >> $outfile
    else
      modfile=`sed 's/\.o\.trs/\.mod/g' <<< $file` >> $outfile
    fi
    echo '|     * '$modfile                        >> $outfile
  done
fi
echo                                               >> $outfile
