// Copyright (C) 2004-2011, Ondra Kamenik

#include "dynare3.h"
#include "dynare_exception.h"
#include "dynare_params.h"

#include "utils/cc/exception.h"
#include "parser/cc/parser_exception.h"
#include "../sylv/cc/SylvException.h"
#include "../kord/random.h"
#include "../kord/global_check.h"
#include "../kord/approximation.h"

int main(int argc, char** argv)
{
	DynareParams params(argc, argv);
	if (params.help) {
		params.printHelp();
		return 0;
	}
	if (params.version) {
		printf("Dynare++ v. %s. Copyright (C) 2004-2011, Ondra Kamenik\n",
			   DYNVERSION);
		printf("Dynare++ comes with ABSOLUTELY NO WARRANTY and is distributed under\n");
		printf("GPL: modules integ, tl, kord, sylv, src, extern and documentation\n");
		printf("LGPL: modules parser, utils\n");
		printf(" for GPL  see http://www.gnu.org/licenses/gpl.html\n");
		printf(" for LGPL see http://www.gnu.org/licenses/lgpl.html\n");
		return 0;
	}
	THREAD_GROUP::max_parallel_threads = params.num_threads;

	try {
		// make journal name and journal
		std::string jname(params.basename);
		jname += ".jnl";
		Journal journal(jname.c_str());

		// make dynare object
		Dynare dynare(params.modname, params.order, params.ss_tol, journal);
		// make list of shocks for which we will do IRFs
        vector<int> irf_list_ind;
		if (params.do_irfs_all)
			for (int i = 0; i < dynare.nexog(); i++)
				irf_list_ind.push_back(i);
		else
			irf_list_ind = ((const DynareNameList&)dynare.getExogNames()).selectIndices(params.irf_list);

		// write matlab files
		FILE* mfd;
		std::string mfile1(params.basename);
		mfile1 += "_f.m";
		if (NULL == (mfd=fopen(mfile1.c_str(), "w"))) {
			fprintf(stderr, "Couldn't open %s for writing.\n", mfile1.c_str());
			exit(1);
		}
		ogdyn::MatlabSSWriter writer0(dynare.getModel(), params.basename.c_str());
		writer0.write_der0(mfd);
		fclose(mfd);

		std::string mfile2(params.basename);
		mfile2 += "_ff.m";
		if (NULL == (mfd=fopen(mfile2.c_str(), "w"))) {
			fprintf(stderr, "Couldn't open %s for writing.\n", mfile2.c_str());
			exit(1);
		}
		ogdyn::MatlabSSWriter writer1(dynare.getModel(), params.basename.c_str());
		writer1.write_der1(mfd);
		fclose(mfd);

		// open mat file
		std::string matfile(params.basename);
		matfile += ".mat";
		mat_t* matfd = Mat_Create(matfile.c_str(), NULL);
		if (matfd == NULL) {
			fprintf(stderr, "Couldn't open %s for writing.\n", matfile.c_str());
			exit(1);
		}

		// write info about the model (dimensions and variables)
		dynare.writeMat(matfd, params.prefix);
		// write the dump file corresponding to the input
		dynare.writeDump(params.basename);


		system_random_generator.initSeed(params.seed);

		tls.init(dynare.order(),
				 dynare.nstat()+2*dynare.npred()+3*dynare.nboth()+
				 2*dynare.nforw()+dynare.nexog());

		Approximation app(dynare, journal, params.num_steps, params.do_centralize, params.qz_criterium);
		try {
			app.walkStochSteady();
		} catch (const KordException& e) {
			// tell about the exception and continue
			printf("Caught (not yet fatal) Kord exception: ");
			e.print();
			JournalRecord rec(journal);
			rec << "Solution routine not finished (" << e.get_message()
				<< "), see what happens" << endrec; 
		}

		std::string ss_matrix_name(params.prefix);
		ss_matrix_name += "_steady_states";
		ConstTwoDMatrix(app.getSS()).writeMat(matfd, ss_matrix_name.c_str());

		// check the approximation
		if (params.check_along_path || params.check_along_shocks
			|| params.check_on_ellipse) {
			GlobalChecker gcheck(app, THREAD_GROUP::max_parallel_threads, journal);
			if (params.check_along_shocks)
				gcheck.checkAlongShocksAndSave(matfd, params.prefix,
											   params.getCheckShockPoints(),
											   params.getCheckShockScale(),
											   params.check_evals);
			if (params.check_on_ellipse)
				gcheck.checkOnEllipseAndSave(matfd, params.prefix,
											 params.getCheckEllipsePoints(),
											 params.getCheckEllipseScale(),
											 params.check_evals);
			if (params.check_along_path)
				gcheck.checkAlongSimulationAndSave(matfd, params.prefix,
												   params.getCheckPathPoints(),
												   params.check_evals);
		}

		// write the folded decision rule to the Mat-4 file
		app.getFoldDecisionRule().writeMat(matfd, params.prefix);

		// simulate conditional
		if (params.num_condper > 0 && params.num_condsim > 0) {
			SimResultsDynamicStats rescond(dynare.numeq(), params.num_condper, 0);
			ConstVector det_ss(app.getSS(),0);
			rescond.simulate(params.num_condsim, app.getFoldDecisionRule(), det_ss, dynare.getVcov(), journal);
			rescond.writeMat(matfd, params.prefix);
		}

		// simulate unconditional
		//const DecisionRule& dr = app.getUnfoldDecisionRule();
		const DecisionRule& dr = app.getFoldDecisionRule();
		if (params.num_per > 0 && params.num_sim > 0) {
			SimResultsStats res(dynare.numeq(), params.num_per, params.num_burn);
			res.simulate(params.num_sim, dr, dynare.getSteady(), dynare.getVcov(), journal);
			res.writeMat(matfd, params.prefix);
			
			// impulse response functions
			if (! irf_list_ind.empty()) {
				IRFResults irf(dynare, dr, res, irf_list_ind, journal);
				irf.writeMat(matfd, params.prefix);
			}
		}

		// simulate with real-time statistics
		if (params.num_rtper > 0 && params.num_rtsim > 0) {
			RTSimResultsStats rtres(dynare.numeq(), params.num_rtper, params.num_burn);
			rtres.simulate(params.num_rtsim, dr, dynare.getSteady(), dynare.getVcov(), journal);
			rtres.writeMat(matfd, params.prefix);
		}

		Mat_Close(matfd);

	} catch (const KordException& e) {
		printf("Caugth Kord exception: ");
		e.print();
		return e.code();
	} catch (const TLException& e) {
		printf("Caugth TL exception: ");
		e.print();
		return 255;
	} catch (SylvException& e) {
		printf("Caught Sylv exception: ");
		e.printMessage();
		return 255;
	} catch (const DynareException& e) {
		printf("Caught Dynare exception: %s\n", e.message());
		return 255;
	} catch (const ogu::Exception& e) {
		printf("Caught ogu::Exception: ");
		e.print();
		return 255;
	} catch (const ogp::ParserException& e) {
		printf("Caught parser exception: %s\n", e.message());
		return 255;
	}

  return 0;
}
