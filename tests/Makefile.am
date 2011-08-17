DYNARE_ROOT = $(abs_top_srcdir)/matlab

MODFILES = \
	ramst.mod \
	ramst_a.mod \
	example1.mod \
	example2.mod \
	example1_use_dll.mod \
	example1_with_tags.mod \
	example1_irf_shocks.mod \
	example1_abs_sign_dynamic.mod \
	t_sgu_ex1.mod \
	osr_example.mod \
	optimal_policy/ramsey.mod \
	optimal_policy/mult_elimination_test.mod \
	discretionary_policy/dennis_1.mod \
	ramst_initval_file.mod \
	ramst_normcdf_and_friends.mod \
	example1_varexo_det.mod \
	predetermined_variables.mod \
	fs2000_nonstationary.mod \
	fs2000_ssfile.mod \
	comments.mod \
	histval_sto.mod \
	histval_det.mod \
	expectation.mod \
	steady_state_operator/standard.mod \
	steady_state_operator/use_dll.mod \
	steady_state_operator/block.mod \
	steady_state_operator/bytecode_test.mod \
	block_bytecode/ireland.mod \
	block_bytecode/ramst_normcdf_and_friends.mod \
	k_order_perturbation/fs2000k2a.mod \
	k_order_perturbation/fs2000k2_use_dll.mod \
	k_order_perturbation/fs2000k_1_use_dll.mod \
	k_order_perturbation/fs2000k3_use_dll.mod \
	k_order_perturbation/fs2000k2_m.mod \
	k_order_perturbation/fs2000k_1_m.mod \
	k_order_perturbation/fs2000k3_m.mod \
	partial_information/PItest3aHc0PCLsimModPiYrVarobsAll.mod \
	partial_information/PItest3aHc0PCLsimModPiYrVarobsCNR.mod \
	arima/mod1.mod \
	arima/mod1a.mod \
	arima/mod1b.mod \
	arima/mod1c.mod \
	arima/mod2.mod \
	arima/mod2a.mod \
	arima/mod2b.mod \
	arima/mod2c.mod \
	fs2000/fs2000.mod \
	fs2000/fs2000a.mod \
	fs2000/fs2000c.mod \
	homotopy/homotopy1_test.mod \
	homotopy/homotopy2_test.mod \
	homotopy/homotopy3_test.mod \
	bvar_a_la_sims/bvar_standalone.mod \
	bvar_a_la_sims/bvar_and_dsge.mod \
	AIM/fs2000x10L9_L.mod \
	AIM/fs2000x10L9_L_AIM.mod \
	AIM/fs2000x10_L9_L.mod \
	AIM/fs2000x10_L9_L_AIM.mod \
	AIM/fs2000_b1L1L.mod \
	AIM/fs2000_b1L1L_AIM.mod \
	AIM/ls2003_2L0L.mod \
	AIM/ls2003_2L0L_AIM.mod \
	AIM/ls2003_2L2L.mod \
	AIM/ls2003_2L2L_AIM.mod \
	conditional_variance_decomposition/example1.mod \
	dsge-var/simul_hybrid.mod \
	dsge-var/dsgevar_forward_calibrated_lambda.mod \
	dsge-var/dsgevar_forward_estimated_lambda.mod \
	external_function/example1_1st_and_2nd_deriv_functions_provided.mod \
	external_function/example1_1st_and_2nd_deriv_functions_provided_dll.mod \
	external_function/example1_1st_deriv_function_provided.mod \
	external_function/example1_1st_deriv_function_provided_dll.mod \
	external_function/example1_no_deriv_functions_provided.mod \
	external_function/example1_no_deriv_functions_provided_dll.mod \
	seeds.mod \
	simul/example1.mod \
	conditional_forecasts/fs2000_cal.mod \
	conditional_forecasts/fs2000_est.mod \
	ms-sbvar/test_exclusions.mod \
	ms-sbvar/test_exclusions_nc.mod \
	ms-sbvar/test_lower_cholesky.mod \
	ms-sbvar/test_lower_cholesky_a.mod \
	ms-sbvar/test_lower_cholesky_nc.mod \
	ms-sbvar/test_upper_cholesky.mod \
	ms-sbvar/test_upper_cholesky_nc.mod \
	ms-sbvar/test1.mod \
	ms-sbvar/test_ms.mod \
	ms-sbvar/test_ms_variances.mod

EXTRA_DIST = \
	$(MODFILES) \
	homotopy/common.mod \
	block_bytecode/ls2003.mod \
	fs2000_ssfile_aux.m \
	printMakeCheckMatlabErrMsg.m \
	printMakeCheckOctaveErrMsg.m \
	ramst_initval_file_data.m \
	run_test.m \
	run_test_matlab.m \
	run_test_octave.m \
	test.m \
	AIM/data_ca1.m \
	AIM/fs2000_b1L1L_AIM_steadystate.m \
	AIM/fs2000_b1L1L_steadystate.m \
	AIM/fsdat.m \
	block_bytecode/run_ls2003.m \
	bvar_a_la_sims/bvar_sample.m \
	external_function/extFunDeriv.m \
	external_function/extFunNoDerivs.m \
	external_function/extFunWithFirstAndSecondDerivs.m \
	fs2000/fs2000a_steadystate.m \
	fs2000/fsdat_simul.m \
	k_order_perturbation/run_fs2000kplusplus.m \
	ls2003/data_ca1.m \
	measurement_errors/data_ca1.m \
	missing/simulate_data_with_missing_observations.m \
	objectives/sgu_ex1.mat \
	conditional_forecasts/fsdat_simul.m

TARGETS = check-test-file

if HAVE_CMD_LINE_MATLAB
TARGETS += check-matlab
endif

if HAVE_OCTAVE
TARGETS += check-octave
endif

check-local: $(TARGETS)

check-test-file: $(MODFILES)
	@set -e; \
	echo "function name=filesToTest()" > filesToTest.m; \
	echo "name = { ..." >> filesToTest.m; \
	for modfile in $(MODFILES); do \
		echo "    '"$$modfile"' ..." >> filesToTest.m; \
	done
	@set -e; \
	echo "    };"  >> filesToTest.m; \
	echo "end"  >> filesToTest.m;

check-matlab:
	$(MATLAB)/bin/matlab -nodesktop -nosplash -r run_test_matlab

check-octave:
	$(OCTAVE) --norc --silent --no-history --path $(DYNARE_ROOT) --eval dynver=\"$(PACKAGE_VERSION)\" run_test_octave.m

clean-local:
	rm -f $(patsubst %.mod, %.m, $(MODS)) \
		$(patsubst %.mod, %_static.*, $(MODS)) \
		$(patsubst %.mod, %_objective_static.m, $(MODS)) \
		$(patsubst %.mod, %_dynamic.*, $(MODS))

	rm -f $(patsubst %.mod, %_results.mat, $(MODS)) \
		$(patsubst %.mod, %_mode.mat, $(MODS)) \
		$(patsubst %.mod, %_mh_mode.mat, $(MODS)) \
		$(patsubst %.mod, %_mean.mat, $(MODS)) \
		$(patsubst %.mod, %_pindx.mat, $(MODS)) \
		$(patsubst %.mod, %_params.mat, $(MODS)) \
		$(patsubst %.mod, %_simul, $(MODS)) \
		$(patsubst %.mod, %.log, $(MODS))

	rm -rf $(patsubst %.mod, %, $(MODS))

	rm -f $(patsubst %.mod, %*.pdf, $(MODS)) \
		$(patsubst %.mod, %*.eps, $(MODS)) \
		$(patsubst %.mod, %*.fig, $(MODS))

	rm -f $(shell find -name g1.mat) \
		$(shell find -name g2.mat) \
		$(shell find -name g3.mat) \
		$(shell find -name H.dat)

	rm -f arima/data1.m arima/data2.m \
		k_order_perturbation/*.jnl \
		k_order_perturbation/*.mat

	rm -rf ramsey_objective

	rm -f fs2000_ssfile_steadystate.m

	rm -f $(shell find -name '*~')

	rm -f dsge-var/datarabanal_hybrid.m

	rm -rf partial_information/PItest3aHc0PCLsimModPiYrVarobsAll_PCL* partial_information/PItest3aHc0PCLsimModPiYrVarobsCNR_PCL*

	rm -rf block_bytecode/ws block_bytecode/ls2003_tmp*
