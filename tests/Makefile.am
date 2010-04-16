DYNARE_ROOT = $(abs_top_srcdir)/matlab

# Under Octave we only test a subset of MOD files, because of missing features (models with unit roots, reading Excel files)
OCTAVE_MODS = \
	ramst.mod \
	ramst_a.mod \
	example1.mod \
	example2.mod \
	example1_use_dll.mod \
	t_sgu_ex1.mod \
	osr_example.mod \
	ramsey.mod \
	ramst_initval_file.mod \
	ramst_normcdf_and_friends.mod \
	ramst_steady_state.mod \
	example1_varexo_det.mod \
	predetermined_variables.mod \
	block_bytecode/fs2000_simk.mod \
	block_bytecode/fs2000_lu.mod \
	block_bytecode/fs2000_bicgstab.mod \
	block_bytecode/fs2000_optpath.mod \
	block_bytecode/fs2000_bytecode.mod \
	block_bytecode/fs2000_block_bytecode.mod \
	block_bytecode/ramst.mod \
	block_bytecode/ireland.mod \
	block_bytecode/ls2003.mod \
	block_bytecode/ramst_normcdf_and_friends.mod \
	block_bytecode/ramst_steady_state.mod \
	k_order_perturbation/fs2000k2a.mod \
	k_order_perturbation/fs2000k2.mod \
	k_order_perturbation/fs2000k_1.mod \
	k_order_perturbation/fs2000k3.mod \
	arima/mod1.mod \
	arima/mod1a.mod \
	arima/mod2.mod \
	arima/mod2c.mod \
	fs2000/fs2000.mod \
	homotopy/homotopy1_test.mod \
	homotopy/homotopy2_test.mod \
	homotopy/homotopy3_test.mod \
	bvar_a_la_sims/bvar_standalone.mod \
	bvar_a_la_sims/bvar_and_dsge.mod \
	AIM/fs2000x10L9_L.mod \
	AIM/fs2000x10L9_L_AIM.mod \
	AIM/fs2000x10_L9_L.mod \
	AIM/fs2000x10_L9_L_AIM.mod \
	conditional_variance_decomposition/example1.mod

MODS = $(OCTAVE_MODS) \
	arima/mod1b.mod \
	arima/mod1c.mod \
	arima/mod2a.mod \
	arima/mod2b.mod \
	fs2000/fs2000a.mod \
	AIM/fs2000_b1L1L.mod \
	AIM/fs2000_b1L1L_AIM.mod \
	AIM/ls2003_2L0L.mod \
	AIM/ls2003_2L0L_AIM.mod \
	AIM/ls2003_2L2L.mod \
	AIM/ls2003_2L2L_AIM.mod \
	block_bytecode/fs2000_gmres.mod \
	block_bytecode/ramst_a.mod

EXTRA_DIST = $(MODS) \
	run_test_octave.m \
	fs2000/fsdat.m \
	test.m \
	objectives \
	ramst_initval_file_data.m \
	homotopy/common.mod \
	block_bytecode/fs2000_common.mod \
	bvar_a_la_sims/bvar_sample.m

TARGETS = check-matlab

if HAVE_OCTAVE
TARGETS += check-octave
endif

check-local: $(TARGETS)

check-octave: $(OCTAVE_MODS)
	@set -e; \
		for modfile in $(OCTAVE_MODS); do \
			$(OCTAVE) --norc --silent --no-history run_test_octave.m $$modfile $(DYNARE_ROOT) $(PACKAGE_VERSION); \
		done
.PHONY: check-octave

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
		k_order_perturbation/*.mat \
		block_bytecode/fs2000_simk.swp

	rm -rf ramsey_objective

	rm -f $(shell find -name '*~')
