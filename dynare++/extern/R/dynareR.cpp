// $Id: dynareR.cpp 862 2006-08-04 17:34:56Z tamas $

// Copyright 2006, Tamas K Papp

#include "dynare3.h"			// Dynare class
#include "approximation.h"		// Approximation class

// exceptions
#include "dynare_exception.h"
#include "parser/cc/parser_exception.h"
#include "utils/cc/exception.h"
#include "SylvException.h"
#include "tl_exception.h"
#include "kord_exception.h"

#include <algorithm>

#include <string.h>

#ifdef DEBUG
#include <stdio.h>
#endif

#include <R_ext/Memory.h>

/** This file containt the C glue functions for an R interface to
 * Dynare++.  Although written in standard C (except for the use of
 * R_alloc), the indexing, calling and memory management conventions
 * of the functions in this file were tailored for R.
 *
 * It is not recommended that you use this interface for anything else
 * but R.
 */ 

/** Error codes: these error codes correspond to possible
 * exceptions. */
#define DYNARER_SYLVEXCEPTION 1
#define DYNARER_DYNAREEXCEPTION 2
#define DYNARER_OGUEXCEPTION 3
#define DYNARER_TLEXCEPTION 4
#define DYNARER_KORDEXCEPTION 5
#define DYNARER_NAMESMATCHINGERROR 6

/** Copies the message into a buffer.  The buffer is allocated and
 * managed by R, ie it will be garbage collected after the .C call
 * returns and the contents are duplicated.
 */
char *passmessage(const char *errormessage) {
	long l = strlen(errormessage);
	char *em = R_alloc(l, 1);
	return strcpy(em, errormessage);
}

/** This function puts the mapping between the newtotal items after
 * nl[offset] and the items in orig into the buffer perm, which has to
 * be at least as long as newtotal.  The function uses R indexing,
 * that is to say, the first index is 1.
 */
int matchnames(const char **orig, int origtotal, 
			   const NameList &nl, int offset, int newtotal,
			   int *perm) {
#ifdef DEBUG
	printf("matching names (R indexing):\n");
#endif
	for (int i=0; i < newtotal; i++) {
		int j;
		for (j=0; j < origtotal; j++)
			if (strcmp(nl.getName(offset+i), *(orig+j))==0) {
				*(perm+i) = j+1;
#ifdef DEBUG
				printf("%d -> %d\n",i+1,j+1);
#endif
				break;
			}
		if (j==origtotal)
			return 1;
	}
	return 0;
}

/** dynareR is the interface function.  The user provides:
 * - a list of endogenous and exogenous variables, a list of
 *   parameters (and the length of each list)
 * - the model equations (modeleq, pointer to a 0-terminated string)
 * - the order of expansion (ord)
 * - journal file name (jnlfile, can be "/dev/null" for no journal)
 * - values for the parametes (parval)
 * - variance-covariance matrix (vcov, stacked by columns, R does
 *   this)
 * - initial values for finding the steady state (initval)
 * - and the number of steps for the approximation algorithm
 *   (numsteps)
 *
 * If successful, the interface will write the results to these
 * buffers:
 * - tensorbuffer for the steady state and the flattened tensors
 * - num_state for the number of endogenous variables that ended up in
 *   the state
 * - mappings to variable names (ordering_state, ordering_endo,
 *   ordering_exo), indices start from 1
 * - the deterministic steady state (newinitval)
 *
 * If dynare throws an exception, the interface tries to catch it and
 * return an error code (error), and error message (errormessage), and
 * if applicable, information on the stability of the model
 * (kordcode).  errormessage is allocated into R's memory, and will be
 * collected after duplication.
 */
extern "C" {
	void dynareR(const char** endo, const int* num_endo,
				 const char** exo, const int* num_exo,
				 const char** par, const int* num_par,
				 const char** equations, const int* ord, const char* jnlfile,
				 const double *parval, const double *vcov, 
				 const double *initval,
				 const int *num_steps,
				 double* tensorbuffer,
				 int *num_state, int *ordering_state,
				 int *ordering_endo, int *ordering_exo,
				 double *newinitval,
				 int* error, char **errormessage, int *kordcode) {
		// construct the model here
		try {	
#ifdef DEBUG					// will print only first var names etc.
			printf("eq: %s\nendo: %d %s\nexo: %d %s\npar: %d %s\nord: %d\n",
				   *equations,*num_endo,*endo,*num_exo,*exo,*num_par,*par,*ord);
#endif
			// create journal
			Journal journal(jnlfile);
			// create Dynare object
			Dynare dynare(endo, *num_endo, exo, *num_exo,
						  par, *num_par, *equations, strlen(*equations),
						  *ord, journal);
			// set Vcov and parameter values
			copy(parval,parval+(*num_par),dynare.getParams().base());
#ifdef DEBUG
			printf("parameter values (%d):\n",dynare.getParams().length());
			dynare.getParams().print();
#endif
			copy(vcov,vcov+(*num_exo)*(*num_exo),dynare.getVcov().base());
#ifdef DEBUG
			printf("vcov matrix:\n");
			dynare.getVcov().print();
#endif
			// set initial values
			Vector iv(initval,*num_endo);
#ifdef DEBUG
			printf("initial values:\n");
			iv.print();
#endif
			dynare.setInitOuter(iv);
			// construct approximation
			tls.init(dynare.order(),
					 dynare.nstat()+2*dynare.npred()+3*dynare.nboth()+
					 2*dynare.nforw()+dynare.nexog());
			Approximation approximation(dynare,journal,*num_steps);
			approximation.walkStochSteady();
			// write the steady state into the buffer
			int ny = dynare.ny();
			const Vector ss(dynare.getSteady());
//			ss = ConstVector(approximation.getSS(), 0); // FIXME allow
//														// for nonzero
			int s = dynare.getStateNames().getNum();
			int sm = s;
			tensorbuffer = copy(ss.base(),ss.base()+ny,tensorbuffer);
			// write the tensors into buffer
			const UnfoldDecisionRule& udr = 
				approximation.getUnfoldDecisionRule();
			for (int i=1; i <= *ord; i++) {
				const UFSTensor* t = udr.get(Symmetry(i));
#ifdef DEBUG
				printf("tensor %d:\n", i);
				t->print();
#endif
				tensorbuffer = copy(t->base(), t->base()+ny*sm, tensorbuffer);
				sm *= s;
				}
			// save number of endogenous states
			*num_state = s-(*num_exo);
			// ordering
#ifdef DEBUG
			printf("all endo names:\n");
			dynare.getAllEndoNames().print();
			printf("all state names:\n");
			dynare.getStateNames().print();
#endif
			if (matchnames(endo, *num_endo, dynare.getAllEndoNames(),
						   0, *num_endo, ordering_endo) ||
				matchnames(endo, *num_endo, dynare.getStateNames(),
						   0, *num_state, ordering_state) ||
				matchnames(exo, *num_exo, dynare.getStateNames(),
						   *num_state, *num_exo, ordering_exo)) {
				*error = DYNARER_NAMESMATCHINGERROR;
				*errormessage = "There was a problem when matching names.  This is weird and should not happen.";
				return;
			}
			// return new init values (first column of SS matrix)
			ConstVector newinit((const GeneralMatrix&) approximation.getSS(), 0);
#ifdef DEBUG
			printf("new initial values:\n");
			newinit.print();
#endif
			copy(newinit.base(),newinit.base()+(*num_endo),newinitval);
		} catch (const SylvException &e) {
			*error = DYNARER_SYLVEXCEPTION;
			char errorbuffer[501];
			e.printMessage(errorbuffer, 500);
			*errormessage = passmessage(errorbuffer);
#ifdef DEBUG
			printf("Caught Sylv exception: ");
			e.printMessage();
#endif
			return;
		} catch (const DynareException &e) {
			*error = DYNARER_DYNAREEXCEPTION;
			*errormessage = passmessage(e.message());
#ifdef DEBUG
			printf("Caught Dynare exception: %s\n", e.message());
#endif
			return;
		}  catch (const ogu::Exception &e) {
			*error = DYNARER_OGUEXCEPTION;
			*errormessage = passmessage(e.message());
#ifdef DEBUG
			printf("Caught ogu::Exception: ");
			e.print();
#endif
			return;
		} catch (const TLException &e) {
			*error = DYNARER_TLEXCEPTION;
			*errormessage = passmessage(e.getmessage());
#ifdef DEBUG
			printf("Caugth TL exception: ");
			e.print();
#endif
			return;
		} catch (const KordException &e) {
			*error = DYNARER_KORDEXCEPTION;
			*errormessage = passmessage(e.getmessage());
			*kordcode = e.code(); // Kord error code
#ifdef DEBUG
			printf("Caugth Kord exception: ");
			e.print();
#endif
			return;
		}
		*error = 0;
		return;}
}
