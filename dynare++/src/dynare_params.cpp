// Copyright (C) 2004-2011, Ondra Kamenik

#include "dynare_params.h"

#include <getopt.h>
#include <cstdio>
#include <cstring>

const char* help_str = 
"usage: dynare++ [--help] [--version] [options] <model file>\n"
"\n"
"    --help               print this message and return\n"
"    --version            print version and return\n"
"\n"
"options:\n"
"    --per <num>          number of periods simulated [100]\n"
"    --sim <num>          number of simulations [80]\n"
"    --rtper <num>        number of RT periods simulated [0]\n"
"    --rtsim <num>        number of RT simulations [0]\n"
"    --condper <num>      number of periods in cond. simulations [0]\n"
"    --condsim <num>      number of conditional simulations [0]\n"
"    --steps <num>        steps towards stoch. SS [0=deter.]\n"
"    --centralize         centralize the rule [do centralize]\n"
"    --no-centralize      do not centralize the rule [do centralize]\n"
"    --prefix <string>    prefix of variables in Mat-4 file [\"dyn\"]\n"
"    --seed <num>         random number generator seed [934098]\n"
"    --order <num>        order of approximation [no default]\n"
"    --threads <num>      number of max parallel threads [2]\n"
"    --ss-tol <num>       steady state calcs tolerance [1.e-13]\n"
"    --check pesPES       check model residuals [no checks]\n"
"                         lower/upper case switches off/on\n"
"                           pP  checking along simulation path\n"
"                           eE  checking on ellipse\n"
"                           sS  checking along shocks\n"
"    --check-evals <num>  max number of evals per residual [1000]\n"
"    --check-num <num>    number of checked points [10]\n"
"    --check-scale <num>  scaling of checked points [2.0]\n"
"    --no-irfs            shuts down IRF simulations [do IRFs]\n"
"    --irfs               performs IRF simulations [do IRFs]\n"
"    --qz-criterium <num> threshold for stable eigenvalues [1.000001]\n"
"\n\n";

// returns the pointer to the first character after the last slash or
// backslash in the string
const char* dyn_basename(const char* str);

DynareParams::DynareParams(int argc, char** argv)
	: modname(NULL), num_per(100), num_sim(80), 
	  num_rtper(0), num_rtsim(0),
	  num_condper(0), num_condsim(0),
	  num_threads(2), num_steps(0),
	  prefix("dyn"), seed(934098), order(-1), ss_tol(1.e-13),
	  check_along_path(false), check_along_shocks(false),
	  check_on_ellipse(false), check_evals(1000), check_num(10), check_scale(2.0),
	  do_irfs_all(true), do_centralize(true), qz_criterium(1.0+1e-6),
	  help(false), version(false)
{
	if (argc == 1 || !strcmp(argv[1],"--help")) {
		help = true;
		return;
	}
	if (argc == 1 || !strcmp(argv[1],"--version")) {
		version = true;
		return;
	}

	modname = argv[argc-1];
	argc--;

	struct option const opts [] = {
		{"periods", required_argument, NULL, opt_per},
		{"per", required_argument, NULL, opt_per},
		{"simulations", required_argument, NULL, opt_sim},
		{"sim", required_argument, NULL, opt_sim},
		{"rtperiods", required_argument, NULL, opt_rtper},
		{"rtper", required_argument, NULL, opt_rtper},
		{"rtsimulations", required_argument, NULL, opt_rtsim},
		{"rtsim", required_argument, NULL, opt_rtsim},
		{"condperiods", required_argument, NULL, opt_condper},
		{"condper", required_argument, NULL, opt_condper},
		{"condsimulations", required_argument, NULL, opt_condsim},
		{"condsim", required_argument, NULL, opt_condsim},
		{"prefix", required_argument, NULL, opt_prefix},
		{"threads", required_argument, NULL, opt_threads},
		{"steps", required_argument, NULL, opt_steps},
		{"seed", required_argument, NULL, opt_seed},
		{"order", required_argument, NULL, opt_order},
		{"ss-tol", required_argument, NULL, opt_ss_tol},
		{"check", required_argument, NULL, opt_check},
		{"check-scale", required_argument, NULL, opt_check_scale},
		{"check-evals", required_argument, NULL, opt_check_evals},
		{"check-num", required_argument, NULL, opt_check_num},
		{"qz-criterium",required_argument, NULL, opt_qz_criterium},
		{"no-irfs", no_argument, NULL, opt_noirfs},
		{"irfs", no_argument, NULL, opt_irfs},
		{"centralize", no_argument, NULL, opt_centralize},
		{"no-centralize", no_argument, NULL, opt_no_centralize},
		{"help", no_argument, NULL, opt_help},
		{"version", no_argument, NULL, opt_version},
		{NULL, 0, NULL, 0}
	};

	int ret;
	int index;
	while (-1 != (ret = getopt_long(argc, argv, "", opts, &index))) {
		switch (ret) {
		case opt_per:
			if (1 != sscanf(optarg, "%d", &num_per))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_sim:
			if (1 != sscanf(optarg, "%d", &num_sim))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_rtper:
			if (1 != sscanf(optarg, "%d", &num_rtper))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_rtsim:
			if (1 != sscanf(optarg, "%d", &num_rtsim))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_condper:
			if (1 != sscanf(optarg, "%d", &num_condper))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_condsim:
			if (1 != sscanf(optarg, "%d", &num_condsim))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_prefix:
			prefix = optarg;
			break;
		case opt_threads:
			if (1 != sscanf(optarg, "%d", &num_threads))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_steps:
			if (1 != sscanf(optarg, "%d", &num_steps))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_seed:
			if (1 != sscanf(optarg, "%d", &seed))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_order:
			if (1 != sscanf(optarg, "%d", &order))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_ss_tol:
			if (1 != sscanf(optarg, "%lf", &ss_tol))
				fprintf(stderr, "Couldn't parse float %s, ignored\n", optarg);
			break;
		case opt_check:
			processCheckFlags(optarg);
			break;
		case opt_check_scale:
			if (1 != sscanf(optarg, "%lf", &check_scale))
				fprintf(stderr, "Couldn't parse float %s, ignored\n", optarg);
			break;
		case opt_check_evals:
			if (1 != sscanf(optarg, "%d", &check_evals))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_check_num:
			if (1 != sscanf(optarg, "%d", &check_num))
				fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
			break;
		case opt_noirfs:
			irf_list.clear();
			do_irfs_all = false;
			break;
		case opt_irfs:
			processIRFList(argc, argv);
			if (irf_list.empty())
				do_irfs_all = true;
			else
				do_irfs_all = false;
			break;
		case opt_centralize:
			do_centralize = true;
			break;
		case opt_no_centralize:
			do_centralize = false;
			break;
		case opt_qz_criterium:
			if (1 != sscanf(optarg, "%lf", &qz_criterium))
				fprintf(stderr, "Couldn't parse float %s, ignored\n", optarg);
			break;
		case opt_help:
			help = true;
			break;
		case opt_version:
			version = true;
			break;
		case '?':
			fprintf(stderr, "Unknown option, ignored\n");
			break;
		}
	}

	// make basename (get rid of the extension)
	basename = dyn_basename(modname);
	std::string::size_type i = basename.rfind('.');
	if (i != std::string::npos)
		basename.erase(i);
}

void DynareParams::printHelp() const
{
	printf("%s", help_str);
}

void DynareParams::processCheckFlags(const char* flags)
{
	for (unsigned int i = 0; i < strlen(flags); i++) {
		switch (flags[i]) {
		case 'p':
			check_along_path = false;
			break;
		case 'P':
			check_along_path = true;
			break;
		case 'e':
			check_on_ellipse = false;
			break;
		case 'E':
			check_on_ellipse = true;
			break;
		case 's':
			check_along_shocks = false;
			break;
		case 'S':
			check_along_shocks = true;
			break;
		default:
			fprintf(stderr, "Unknown check type selection character <%c>, ignored.\n", flags[i]);
		}
	}
}

void DynareParams::processIRFList(int argc, char** argv)
{
	irf_list.clear();
	while (optind < argc && *(argv[optind]) != '-') {
		irf_list.push_back(argv[optind]);
		optind++;
	}
}

const char* dyn_basename(const char* str)
{
	int i = strlen(str);
	while (i > 0 && str[i-1] != '/' && str[i-1] != '\\')
		i--;
	return str+i;
}

// Local Variables:
// mode:C++
// End:
