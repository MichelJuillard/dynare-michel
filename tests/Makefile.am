DYNARE_ROOT = $(abs_top_srcdir)/matlab

MODS = \
	ramst.mod \
	ramst_a.mod \
	example1.mod \
	example2.mod \
	example1_use_dll.mod \
	example1_with_tags.mod\
	t_sgu_ex1.mod \
	osr_example.mod \
	optimal_policy/ramsey.mod \
	optimal_policy/mult_elimination_test.mod \
	ramst_initval_file.mod \
	ramst_normcdf_and_friends.mod \
	example1_varexo_det.mod \
	predetermined_variables.mod \
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
	k_order_perturbation/fs2000k2.mod \
	k_order_perturbation/fs2000k_1.mod \
	k_order_perturbation/fs2000k3.mod \
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
	seeds.mod

EXTRA_DIST = $(MODS) \
	run_test_octave.m \
	block_bytecode/run_block_bytecode_tests.m \
	block_bytecode/run_ls2003.m \
	block_bytecode/ls2003.mod \
	fs2000/fsdat_simul.m \
	test.m \
	objectives \
	ramst_initval_file_data.m \
	homotopy/common.mod \
	bvar_a_la_sims/bvar_sample.m

TARGETS = check-matlab

if HAVE_OCTAVE
TARGETS += check-octave check-block-bytecode
endif

check-local: $(TARGETS)

check-octave: $(MODS)
	@set -e; \
		for modfile in $(MODS); do \
			$(OCTAVE) --norc --silent --no-history run_test_octave.m $$modfile $(DYNARE_ROOT) $(PACKAGE_VERSION); \
		done
.PHONY: check-octave

check-block-bytecode:
	cd block_bytecode && $(OCTAVE) --norc --silent --no-history run_block_bytecode_tests.m $(DYNARE_ROOT) $(PACKAGE_VERSION)
.PHONY: check-block-bytecode

check-matlab: $(MODS)
# MATLAB stuff to be added here
.PHONY: check-matlab

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

	rm -f dsge-var/datarabanal_hybrid.mat

	rm -rf partial_information/PItest3aHc0PCLsimModPiYrVarobsAll_PCL* partial_information/PItest3aHc0PCLsimModPiYrVarobsCNR_PCL*

	rm -rf block_bytecode/ws block_bytecode/ls2003_tmp*
