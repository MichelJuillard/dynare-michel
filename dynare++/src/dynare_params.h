// $Id: dynare_params.h 2347 2009-03-24 11:54:29Z kamenik $

// Copyright 2004, Ondra Kamenik

/*
along shocks: m    mult    max_evals
ellipse:      m    mult    max_evals  (10*m) (0.5*mult)
simul:        m            max_evals  (10*m)

--check-scale 2.0 --check-evals 1000 --check-num 10 --check PES
 */

#include <vector>
#include <string>

struct DynareParams {
	const char* modname;
	std::string basename;
	int num_per;
	int num_burn;
	int num_sim;
	int num_rtper;
	int num_rtsim;
	int num_condper;
	int num_condsim;
	int num_threads;
	int num_steps;
	const char* prefix;
	int seed;
	int order;
	/** Tolerance used for steady state calcs. */
	double ss_tol;
	bool check_along_path;
	bool check_along_shocks;
	bool check_on_ellipse;
	int check_evals;
	int check_num;
	double check_scale;
	/** Flag for doing IRFs even if the irf_list is empty. */
	bool do_irfs_all;
	/** List of shocks for which IRF will be calculated. */
	std::vector<const char*> irf_list;
	bool do_centralize;
	double qz_criterium;
	bool help;
	bool version;
	DynareParams(int argc, char** argv);
	void printHelp() const;
	int getCheckShockPoints() const
		{return check_num;}
	double getCheckShockScale() const
		{return check_scale;}
	int getCheckEllipsePoints() const
		{return 10*check_num;}
	double getCheckEllipseScale() const
		{return 0.5*check_scale;}
	int getCheckPathPoints() const
		{return 10*check_num;}
private:
	enum {opt_per, opt_burn, opt_sim, opt_rtper, opt_rtsim, opt_condper, opt_condsim,
		  opt_prefix, opt_threads,
		  opt_steps, opt_seed, opt_order, opt_ss_tol, opt_check,
		  opt_check_along_path, opt_check_along_shocks, opt_check_on_ellipse,
		  opt_check_evals, opt_check_scale, opt_check_num, opt_noirfs, opt_irfs,
                  opt_help, opt_version, opt_centralize, opt_no_centralize, opt_qz_criterium};
	void processCheckFlags(const char* flags);
	/** This gathers strings from argv[optind] and on not starting
	 * with '-' to the irf_list. It stops one item before the end,
	 * since this is the model file. */  
	void processIRFList(int argc, char** argv);
};

// Local Variables:
// mode:C++
// End:
