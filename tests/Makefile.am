.NOTPARALLEL: check-local

DYNARE_ROOT = $(abs_top_srcdir)/matlab

TESTS_ENVIRONMENT = \
	$(OCTAVE) --norc --silent --no-history --path $(DYNARE_ROOT) --eval dynver=\"$(PACKAGE_VERSION)\"

TESTS = \
	run_test_octave.m

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
