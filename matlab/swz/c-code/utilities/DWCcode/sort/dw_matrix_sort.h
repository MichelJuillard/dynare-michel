
#ifndef __SORT_MATRICES__
#define __SORT_MATRICES__

#include "swzmatrix.h"

TVector SortVectorAscending(TVector x, TVector y);
TVector SortVectorDescending(TVector x, TVector y);
TMatrix SortMatrixRowsAscending(TMatrix X, TMatrix Y, int j);
TMatrix SortMatrixRowsDescending(TMatrix X, TMatrix Y, int j);
TMatrix SortMatrixColumnsAscending(TMatrix X, TMatrix Y, int i);
TMatrix SortMatrixColumnsDescending(TMatrix X, TMatrix Y, int i);

#endif
