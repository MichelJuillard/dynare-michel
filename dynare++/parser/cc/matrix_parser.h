// Copyright (C) 2006, Ondra Kamenik

// $Id: matrix_parser.h 762 2006-05-22 13:00:07Z kamenik $

#ifndef OGP_MATRIX_PARSER
#define OGP_MATRIX_PARSER

#include <cstdlib> // For NULL
#include <vector>

namespace ogp {
	using std::vector;

	/** This class reads the given string and parses it as a
	 * matrix. The matrix is read row by row. The row delimiter is
	 * either a newline character or semicolon (first newline
	 * character after the semicolon is ignored), the column delimiter
	 * is either blank character or comma. A different number of items
	 * in the row is not reconciliated, we do not construct a matrix
	 * here. The class provides only an iterator to go through all
	 * read items, the iterator provides information on row number and
	 * column number of the item. */
	class MPIterator;
	class MatrixParser {
		friend class MPIterator;
	protected:
		/** Raw data as they were read. */
		vector<double> data;
		/** Number of items in each row. */
		vector<int> row_lengths;
		/** Maximum number of row lengths. */
		int nc;
	public:
		MatrixParser()
			: nc(0) {}
		MatrixParser(const MatrixParser& mp)
			: data(mp.data), row_lengths(mp.row_lengths), nc(mp.nc) {}
		virtual ~MatrixParser() {}
		/** Return a number of read rows. */
		int nrows() const
			{return (int) row_lengths.size();}
		/** Return a maximum number of items in the rows. */
		int ncols() const
			{return nc;}
		/** Parses a given data. This initializes the object data. */
		void parse(int length, const char* stream);
		/** Adds newly read item. This should be called from bison
		 * parser. */
		void add_item(double v);
		/** Starts a new row. This should be called from bison
		 * parser. */
		void start_row();
		/** Process a parse error from the parser. */
		void error(const char* mes) const;
		/** Return begin iterator. */
		MPIterator begin() const;
		/** Return end iterator. */
		MPIterator end() const;
	protected:
		/** Returns an index of the first non-empty row starting at
		 * start. If the start row is non-empty, returns the start. If
		 * there is no other non-empty row, returns
		 * row_lengths.size(). */
		int find_first_non_empty_row(int start = 0) const;
	};

	/** This is an iterator intended to iterate through a matrix parsed
	 * by MatrixParser. The iterator provides only read-only access. */
	class MPIterator {
		friend class MatrixParser;
	protected:
		/** Reference to the matrix parser. */
		const MatrixParser* p;
		/** The index of the pointed item in the matrix parser. */
		unsigned int i;
		/** The column number of the pointed item starting from zero. */
		int c;
		/** The row number of the pointed item starting from zero. */
		int r;

	public:
		MPIterator() : p(NULL), i(0), c(0), r(0) {}
		/** Constructs an iterator pointing to the beginning of the
		 * parsed matrix. */
		MPIterator(const MatrixParser& mp);
		/** Constructs an iterator pointing to the past-the-end of the
		 * parsed matrix. */
		MPIterator(const MatrixParser& mp, const char* dummy);
		/** Return read-only reference to the pointed item. */
		const double& operator*() const
			{return p->data[i];}
		/** Return a row index of the pointed item. */
		int row() const
			{return r;}
		/** Return a column index of the pointed item. */
		int col() const
			{return c;}
		/** Assignment operator. */
		const MPIterator& operator=(const MPIterator& it)
			{p = it.p; i = it.i; c = it.c; r = it.r; return *this;}
		/** Return true if the iterators are the same, this is if they
		 * have the same underlying object and the same item index. */ 
		bool operator==(const MPIterator& it) const
			{return it.p == p && it.i == i;}
		/** Negative of the operator==. */
		bool operator!=(const MPIterator& it) const
			{return ! (it == *this);} 
		/** Increment operator. */
		MPIterator& operator++();
	};
};



#endif

// Local Variables:
// mode:C++
// End:
