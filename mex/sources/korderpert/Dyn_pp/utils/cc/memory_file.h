// Copyright (C) 2005, Ondra Kamenik

// $Id: memory_file.h 762 2006-05-22 13:00:07Z kamenik $

#ifndef OGU_MEMORY_FILE
#define OGU_MEMORY_FILE

namespace ogu {
	/** This function calculates an offset of a given position in a
	 * given string. The position is given by the line number and by
	 * the offset in the line (both starting from 1). */
	int calc_pos_offset(int length, const char* str, int line, int col);
	/** This function calculates a line number and column number of a
	 * character given by the offset in the string. It is inverse to
	 * calc_pos_offset. */
	void calc_pos_line_and_col(int length, const char* str, int offset,
							   int& line, int& col);

	/** This class opens a given file and makes its copy in memory and
	 * appends it with the '\0' character. Since the type of length is
	 * int, it can store files with size at most 4GB. If the file
	 * could be opened for reading, data is NULL and length is -1. If
	 * the file is empty but exists, len is zero and data points to a
	 * newly allocated memory containing '\0' character at the end. */
	class MemoryFile {
	protected:
		int len;
		char* data;
	public:
		MemoryFile(const char* fname);
		virtual ~MemoryFile()
			{if (data) delete [] data;}
		int length() const
			{return len;}
		const char* base() const
			{return data;}
		bool exists() const
			{return len != -1;}
		/** Return the offset of a character in the given line
		 * (starting from 1) with the given offset in the line. */
		int offset(int line, int lineoff) const
			{return calc_pos_offset(len, data, line, lineoff);}
		/** Return the line number and column number of the character
		 * defined by the offset. */
		void line_and_col(int offset, int& line, int& col) const
			{calc_pos_line_and_col(len, data, offset, line, col);}
	};

	
};


#endif

// Local Variables:
// mode:C++
// End:
