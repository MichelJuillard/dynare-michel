#include "pascal_triangle.h"
#include <cstdio>

using namespace ogu;

PascalTriangle ptriang;

void PascalRow::setFromPrevious(const PascalRow& prev)
{
	k = prev.k + 1;
	clear();
	prolong(prev);
}

/** This prolongs the PascalRow. If it is empty, we set the first item
 * to k+1, which is noverk(k+1,k) which is the second item in the real
 * pascal row, which starts from noverk(k,k)=1. Then we calculate
 * other items from the provided row which must be the one with k-1.*/
void PascalRow::prolong(const PascalRow& prev)
{
	if (size() == 0)
		push_back(k+1);
	int last = back();
	for (unsigned int i = size(); i < prev.size(); i++) {
		last += prev[i];
		push_back(last);
	}
}

void PascalRow::prolongFirst(int n)
{
	// todo: check n = 1;
	for (int i = (int)size()+2; i <= n; i++)
		push_back(i);
}

void PascalRow::print() const
{
	printf("k=%d\n",k);
	for (unsigned int i = 0; i < size(); i++)
		printf("%d ",operator[](i));
	printf("\n");
}

int PascalTriangle::max_n() const
{
	return (int)(tr[0].size()+1);
}

int PascalTriangle::max_k() const
{
	return (int)tr.size();
}

void PascalTriangle::ensure(int n, int k)
{
	// add along n
	if (n > max_n()) {
		tr[0].prolongFirst(n);
		for (int i = 2; i <= max_k(); i++)
			tr[i-1].prolong(tr[i-2]);
	}

	if (k > max_k()) {
		for (int i = max_k()+1; i <= k; i++) {
			PascalRow r;
			tr.push_back(r);
			tr.back().setFromPrevious(tr[i-2]);
		}
	}
}

int PascalTriangle::noverk(int n, int k)
{
	// todo: rais if out of bounds
	if (n-k < k)
		k = n-k;
	if (k == 0)
		return 1;
	ensure(n, k);
	return (tr[k-1])[n-1-k];
}

void PascalTriangle::print() const
{
	for (unsigned int i = 0; i < tr.size(); i++)
		tr[i].print();
}
