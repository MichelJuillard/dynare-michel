// Copyright (C) 2005, Ondra Kamenik

// $Id: pascal_triangle.h 762 2006-05-22 13:00:07Z kamenik $

#ifndef PASCAL_TRIANGLE_H
#define PASCAL_TRIANGLE_H

#include <vector>

namespace ogu {

	using std::vector;

	class PascalRow : public vector<int> {
		int k;
	public:
		PascalRow()
			: vector<int>(), k(1)
			{ push_back(2); }
		void setFromPrevious(const PascalRow& prev);
		void prolong(const PascalRow& prev);
		void prolongFirst(int n);
		void print() const;
	};

	class PascalTriangle {
		vector<PascalRow> tr;
	public:
		PascalTriangle()
			{tr.push_back(PascalRow());}
		PascalTriangle(const PascalTriangle& triang)
			: tr(triang.tr) {}
		const PascalTriangle& operator=(const PascalTriangle& triang)
			{ tr = triang.tr; return *this;}
		int noverk(int n, int k);
		void print() const;
	protected:
		void ensure(int n, int k);
		int max_n() const;
		int max_k() const;
	};
};

extern ogu::PascalTriangle ptriang;


#endif

// Local Variables:
// mode:C++
// End:
