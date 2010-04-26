
#include "swzmatrix.h"

char* SphericalType(void);
void SetupSpherical_Gaussian(int n);
void SetupSpherical_Uniform(int n);
void SetupSpherical_Power(int n, PRECISION k);
void SetupSpherical_TruncatedPower(int n, PRECISION k, PRECISION a);
void SetupSpherical_Table(int n, PRECISION *table, int m);
void SetupSpherical_TruncatedGaussian(int n, PRECISION r1, PRECISION r2);

PRECISION DrawSpherical(TVector x);
PRECISION LogSphericalDensity(PRECISION r);
TVector SphericalCumulativeDensity(PRECISION max, int bins, int cum_bins);

void TestSpherical(FILE *f, char *filename, PRECISION max);

PRECISION UniformUnitBall(TVector x);
PRECISION PowerUnitBall(TVector x, PRECISION k);
PRECISION TruncatedPowerUnitBall(TVector x, PRECISION k, PRECISION a);
PRECISION SphericalTable(TVector x, PRECISION *table, int m);
PRECISION LogSphericalTableDensity(PRECISION r, PRECISION *table, int m);

