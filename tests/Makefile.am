DYNARE_ROOT = $(abs_top_srcdir)/matlab

# Under Octave we only test a subset of MOD files, because of missing features (models with unit roots, reading Excel files)
OCTAVE_MODS = \
	ramst.mod \
	ramst_a.mod \
	example1.mod \
	t_sgu_ex1.mod \
	ramsey.mod \
	ramst_initval_file.mod \
	example1_varexo_det.mod \
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
	AIM/fs2000x10_L9_L_AIM.mod

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
	AIM/ls2003_2L2L_AIM.mod

EXTRA_DIST = $(MODS) \
	run_test_octave.m \
	fs2000/fsdat.m \
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
		$(patsubst %.mod, %_static.m, $(MODS)) \
		$(patsubst %.mod, %_objective_static.m, $(MODS)) \
		$(patsubst %.mod, %_dynamic.m, $(MODS))

	rm -f $(patsubst %.mod, %_results.mat, $(MODS)) \
		$(patsubst %.mod, %_mode.mat, $(MODS)) \
		$(patsubst %.mod, %_mh_mode.mat, $(MODS)) \
		$(patsubst %.mod, %_mean.mat, $(MODS)) \
		$(patsubst %.mod, %_pindx.mat, $(MODS)) \
		$(patsubst %.mod, %_params.mat, $(MODS)) \
		$(patsubst %.mod, %.log, $(MODS))
	rm -rf $(patsubst %.mod, %, $(MODS))

	rm -f $(patsubst %.mod, %*.pdf, $(MODS)) \
		$(patsubst %.mod, %*.eps, $(MODS)) \
		$(patsubst %.mod, %*.fig, $(MODS))

	rm -f $(shell find -name g1.mat) \
		$(shell find -name g2.mat) \
		$(shell find -name g3.mat) \
		$(shell find -name H.dat)

	rm -f arima/data1.m arima/data2.m

	rm -f $(shell find -name '*~')
