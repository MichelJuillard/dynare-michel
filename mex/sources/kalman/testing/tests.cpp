// $Id: tests.cpp 534 2005-11-30 13:58:11Z kamenik $
// Copyright 2005, Ondra Kamenik

#include "../cc/kalman.h"
#include "../cc/ts_exception.h"
#include "ascii_matrix.h"

#include "GeneralMatrix.h"
#include "Vector.h"
#include "SylvException.h"

#include <sys/time.h>
#include <math.h>


// gettimeofday for MinGW
#ifdef __MINGW32__
#define _W32_FT_OFFSET (116444736000000000LL) 

typedef struct _filetime {
	unsigned long dwLowDateTime;
	unsigned long dwHighDateTime;
} filetime;

extern "C" {
	void __stdcall GetSystemTimeAsFileTime(filetime*);
};

typedef union {
	long long ns100; // time since 1 Jan 1601 in 100ns units
	filetime ft;
} w32_ftv;

void D_gettimeofday(struct timeval* p, struct timezone* tz)
{
	w32_ftv _now;
	GetSystemTimeAsFileTime( &(_now.ft) );
	p->tv_usec=(long)((_now.ns100 / 10LL) % 1000000LL );
	p->tv_sec= (long)((_now.ns100-_W32_FT_OFFSET)/10000000LL);
	return;
}

#else
#define D_gettimeofday gettimeofday
#endif // gettimeofday for MinGW


struct AsciiKalmanTask {
	AsciiMatrix Z;
	AsciiMatrix H;
	AsciiMatrix T;
	AsciiMatrix R;
	AsciiMatrix Q;
	AsciiMatrix Pstar;
	AsciiMatrix Pinf;
	AsciiMatrix a;
	AsciiMatrix Y;
	AsciiKalmanTask(const char* prefix)
		: Z(std::string(prefix) + "_Z.dat"),
		  H(std::string(prefix) + "_H.dat"),
		  T(std::string(prefix) + "_T.dat"),
		  R(std::string(prefix) + "_R.dat"),
		  Q(std::string(prefix) + "_Q.dat"),
		  Pstar(std::string(prefix) + "_Pstar.dat"),
		  Pinf(std::string(prefix) + "_Pinf.dat"),
		  a(std::string(prefix) + "_a.dat"),
		  Y(std::string(prefix) + "_Y.dat")
		{}
};

// WallTimer class. Constructor saves the wall time, destructor
// cancels the current time from the saved, and prints the message
// with time information
class WallTimer {
	char mes[100];
	struct timeval start;
	bool new_line;
public:
	WallTimer(const char* m, bool nl = true)
		{strcpy(mes, m);new_line = nl; D_gettimeofday(&start, NULL);}
	~WallTimer()
		{
			struct timeval end;
			D_gettimeofday(&end, NULL);
			printf("%s%8.4g", mes,
				   end.tv_sec-start.tv_sec + (end.tv_usec-start.tv_usec)*1.0e-6);
			if (new_line)
				printf("\n");
		}
};

/****************************************************/
/*     declaration of TestRunnable class            */
/****************************************************/
class TestRunnable {
	char name[100];
public:
	TestRunnable(const char* n)
		{strncpy(name, n, 100);}
	bool test() const;
	virtual bool run() const =0;
	const char* getName() const
		{return name;}
protected:
	static bool filter_and_smoother(const char* prefix, bool diffuse_flag);
};

bool TestRunnable::test() const
{
	printf("Running test <%s>\n",name);
	bool passed;
	{
		WallTimer tim("Wall clock time ", false);
		passed = run();
	}
	if (passed) {
		printf("............................ passed\n\n");
		return passed;
	} else {
		printf("............................ FAILED\n\n");
		return passed;
	}
}

/****************************************************/
/*     definition of TestRunnable static methods    */
/****************************************************/
bool TestRunnable::filter_and_smoother(const char* prefix, bool diffuse_flag)
{
	AsciiKalmanTask akt(prefix);
	StateInit* init;
	if (diffuse_flag)
		init = new StateInit(akt.Pstar, akt.Pinf, akt.a.getData());
	else
		init = new StateInit(akt.Pstar, akt.a.getData());

	KalmanTask kt(akt.Y, akt.Z, akt.H, akt.T, akt.R, akt.Q, *init);

	// multivariate
	int per;
	int d;
	double ll;
	GeneralMatrix alpha(akt.T.numRows(), akt.Y.numCols());
	GeneralMatrix eta(akt.R.numCols(), akt.Y.numCols());
	GeneralMatrix V(akt.T.numRows(), akt.T.numRows()*akt.Y.numCols());
	SmootherResults sres(akt.Y.numCols());
	{
		WallTimer tim("\tMultivariate time      ", true);
		ll = kt.filter_and_smooth(sres, per, d);
		printf("\t\tll=%f per=%d d=%d\n", ll, per, d);
		if (per == akt.Y.numCols()) {
			sres.exportAlpha(alpha);
			sres.exportEta(eta);
			sres.exportV(V);
		} else {
			printf("\t\tNot finished.\n");
		}
	}

	// univariate
	KalmanUniTask kut(kt);
	int per1;
	int d1;
	double ll1;
	GeneralMatrix alpha1(akt.T.numRows(), akt.Y.numCols());
	GeneralMatrix eta1(akt.R.numCols(), akt.Y.numCols());
	GeneralMatrix V1(akt.T.numRows(), akt.T.numRows()*akt.Y.numCols());
	SmootherResults sres1(akt.Y.numCols()*akt.Y.numRows());
	{
		WallTimer tim("\tUnivariate time        ", true);
		int dd;
		ll1 = kut.filter_and_smooth(sres1, per1, dd);
		per1 /= akt.Y.numRows();
		d1 = dd/akt.Y.numRows();
		printf("\t\tll=%f per=%d d=%d(%d)\n", ll1, per1, d1, dd);
		if (per1 == akt.Y.numCols()) {
			SmootherResults sres_uni(akt.Y.numCols());
			sres_uni.import(sres1, akt.Y.numRows());
			sres_uni.exportAlpha(alpha1);
			sres_uni.exportEta(eta1);
			sres_uni.exportV(V1);
		} else {
			printf("\t\tNot finished.\n");
		}
	}

	// compare
	if (per == per1 && per == akt.Y.numCols()) {
		WallTimer tim("\tComparison time        ", true);
		alpha.add(-1.0, alpha1);
		eta.add(-1.0, eta1);
		V.add(-1.0, V1);
		int maxd = std::max(d,d1);
		for (int t = 1; t <= maxd; t++) {
			Vector alphat(alpha, t-1);
			printf("\t\tt=%d   alpha error %10.6g\n",t,alphat.getMax());
			Vector etat(eta, t-1);
			printf("\t\tt=%d   eta   error %10.6g\n",t,etat.getMax());
			GeneralMatrix Vt(V, 0, (t-1)*akt.T.numRows(), akt.T.numRows(), akt.T.numRows());
			printf("\t\tt=%d   V     error %10.6g\n",t,V.getData().getMax());
		}
		GeneralMatrix alpha_rest(alpha, 0, maxd, akt.T.numRows(), alpha.numCols()-maxd);
		printf("\t\tt=%d.. alpha error %10.6g\n",maxd+1,alpha_rest.getData().getMax());
		GeneralMatrix eta_rest(eta, 0, maxd, akt.R.numCols(), eta.numCols()-maxd);
		printf("\t\tt=%d.. eta   error %10.6g\n",maxd+1,eta_rest.getData().getMax());
		GeneralMatrix V_rest(V, 0, maxd*akt.T.numRows(), akt.T.numRows(),
							 V.numCols()-maxd*akt.T.numRows());
		printf("\t\tt=%d.. V     error %10.6g\n",maxd+1,V_rest.getData().getMax());
	}

	delete init;

	return true;
}


/****************************************************/
/*     definition of TestRunnable subclasses        */
/****************************************************/
class SmallNonDiffuse : public TestRunnable {
public:
	SmallNonDiffuse()
		: TestRunnable("Non-diffuse small (p=2,m=3,r=4)") {}

	bool run() const
		{
			filter_and_smoother("small2x3x4", false);
			return true;
		}
};

class SmallDiffuse : public TestRunnable {
public:
	SmallDiffuse()
		: TestRunnable("Diffuse small (p=2,m=3,r=4)") {}

	bool run() const
		{
			return filter_and_smoother("small2x3x4", true);
		}
};

class MiddleNonDiffuse : public TestRunnable {
public:
	MiddleNonDiffuse()
		: TestRunnable("Non-diffuse middle (p=10,m=15,r=12)") {}

	bool run() const
		{
			return filter_and_smoother("10x15x12", false);
		}
};

class MiddleDiffuse : public TestRunnable {
public:
	MiddleDiffuse()
		: TestRunnable("Diffuse middle (p=10,m=15,r=12)") {}

	bool run() const
		{
			return filter_and_smoother("10x15x12", true);
		}
};

class SOEDiffuse : public TestRunnable {
public:
	SOEDiffuse()
		: TestRunnable("Diffuse soe (p=8,m=25,r=15)") {}

	bool run() const
		{
			return filter_and_smoother("soe8x25x15", true);
		}
};

int main()
{
	TestRunnable* all_tests[50];
	// fill in vector of all tests
	int num_tests = 0;
	all_tests[num_tests++] = new SmallNonDiffuse();
	all_tests[num_tests++] = new SmallDiffuse();
	all_tests[num_tests++] = new MiddleNonDiffuse();
	all_tests[num_tests++] = new MiddleDiffuse();
	all_tests[num_tests++] = new SOEDiffuse();

	// launch the tests
	int success = 0;
	for (int i = 0; i < num_tests; i++) {
		try {
			if (all_tests[i]->test())
				success++;
		} catch (const TSException& e) {
			printf("Caugth TS exception in <%s>:\n", all_tests[i]->getName());
			e.print();
		} catch (SylvException& e) {
			printf("Caught Sylv exception in <%s>:\n", all_tests[i]->getName());
			e.printMessage();
		}
	}

	printf("There were %d tests that failed out of %d tests run.\n",
		   num_tests - success, num_tests);

	// destroy
	for (int i = 0; i < num_tests; i++) {
		delete all_tests[i];
	}

	return 0;
}
