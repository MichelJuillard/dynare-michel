// Copyright (C) 2006, Ondra Kamenik

// $Id: atom_assignings.cpp 92 2007-04-19 11:38:21Z ondra $

#include "atom_assignings.h"
#include "location.h"
#include "parser_exception.h"

#include "utils/cc/exception.h"

#include <limits>
#include <iostream>

using namespace ogp;

AtomAssignings::AtomAssignings(const AtomAssignings& aa, ogp::StaticAtoms& a)
	: atoms(a), expr(aa.expr, atoms), left_names(aa.left_names),
	  order(aa.order)
{
	// fill the lname2expr
	for (Tvarintmap::const_iterator it = aa.lname2expr.begin();
		 it != aa.lname2expr.end(); ++it)
		lname2expr.insert(Tvarintmap::value_type(left_names.query((*it).first), (*it).second));
}

/** A global symbol for passing info to the AtomAssignings from
 * asgn_parse(). */
AtomAssignings* aparser;

/** The declaration of functions defined in asgn_ll.cc and asgn_tab.cc
 * generated from assign.lex assign.y */
void* asgn__scan_buffer(char*, size_t);
void asgn__destroy_buffer(void*);
void asgn_parse();
extern location_type asgn_lloc;

void AtomAssignings::parse(int length, const char* stream)
{
	char* buffer = new char[length+2];
	strncpy(buffer, stream, length);
	buffer[length] = '\0';
	buffer[length+1] = '\0';
	asgn_lloc.off = 0;
	asgn_lloc.ll = 0;
	void* p = asgn__scan_buffer(buffer, (unsigned int)length+2);
	aparser = this;
	asgn_parse();
	delete [] buffer;
	asgn__destroy_buffer(p);
}

void AtomAssignings::error(const char* mes)
{
	throw ParserException(mes, asgn_lloc.off);
}

void AtomAssignings::add_assignment_to_double(const char* name, double val)
{
	// if left hand side is a registered atom, insert it to tree
	int t;
	try {
		if (atoms.check(name))
			t = expr.add_nulary(name);
		else
			t = -1;
	} catch (const ParserException& e) {
		t = -1;
	}
	// register left hand side in order
	order.push_back(t);

	// add the double to the tree
	char tmp[100];
	sprintf(tmp, "%30.25g", val);
	try {
		expr.parse(strlen(tmp), tmp);
	} catch (const ParserException& e) {
		// should never happen
		throw ParserException(string("Error parsing double ")+tmp+": "+e.message(), 0);
	}

	// register name of the left hand side and put to lname2expr
	const char* ss = left_names.insert(name);
	lname2expr.insert(Tvarintmap::value_type(ss, order.size()-1));
}

void AtomAssignings::add_assignment(int asgn_off, const char* str, int name_len,
									int right_off, int right_len)
{
	// the order of doing things here is important: since the
	// FormulaParser requires that all references from the i-th tree
	// refere to trees with index lass than i, so to capture also a
	// nulary term for the left hand side, it must be inserted to the
	// expression tree before the expression is parsed.

	// find the name in the atoms, make copy of name to be able to put
	// '\0' at the end
	char* buf = new char[name_len+1];
	strncpy(buf, str, name_len);
	buf[name_len] = '\0';
	// if left hand side is a registered atom, insert it to tree
	int t;
	try {
		t = atoms.check(buf);
		if (t == -1)
			t = expr.add_nulary(buf);
	} catch (const ParserException& e) {
		atoms.register_name(buf);
		t = expr.add_nulary(buf);
	}
	// register left hand side in order
	order.push_back(t);

	// parse expression on the right
	try {
		expr.parse(right_len, str+right_off);
	} catch (const ParserException& e) {
		throw ParserException(e, asgn_off+right_off);
	}

	// register name of the left hand side and put to lname2expr
	const char* ss = left_names.insert(buf);
	if (lname2expr.find(ss) != lname2expr.end()) {
		// Prevent the occurrence of #415
		std::cerr << "Changing the value of " << ss << " is not supported. Aborting." << std::endl;
		exit(EXIT_FAILURE);
	}
	lname2expr[ss] = order.size()-1;

	// delete name
	delete [] buf;
}

void AtomAssignings::apply_subst(const AtomSubstitutions::Toldnamemap& mm)
{
	// go through all old variables and see what are their derived new
	// variables
	for (AtomSubstitutions::Toldnamemap::const_iterator it = mm.begin();
		 it != mm.end(); ++it) {
		const char* oldname = (*it).first;
		const AtomSubstitutions::Tshiftnameset& sset = (*it).second;
		if (! sset.empty()) {
			int told = atoms.index(oldname);
			if (told < 0 && ! atoms.get_name_storage().query(oldname))
				atoms.register_name(oldname);
			if (told == -1)
				told = expr.add_nulary(oldname);
			// at least one substitution here, so make an expression
			expr.add_formula(told);
			// say that this expression is not assigned to any atom
			order.push_back(-1);
			// now go through all new names derived from the old name and
			// reference to the newly added formula
			for (AtomSubstitutions::Tshiftnameset::const_iterator itt = sset.begin();
				 itt != sset.end(); ++itt) {
				const char* newname = (*itt).first;
				const char* nn = left_names.insert(newname);
				lname2expr.insert(Tvarintmap::value_type(nn, expr.nformulas()-1));
			}
		}
	}
}

void AtomAssignings::print() const
{
	printf("Atom Assignings\nExpressions:\n");
	expr.print();
	printf("Left names:\n");
	for (Tvarintmap::const_iterator it = lname2expr.begin();
		 it != lname2expr.end(); ++it)
		printf("%s ==> %d (t=%d)\n", (*it).first, expr.formula((*it).second), order[(*it).second]);
}

void AtomAsgnEvaluator::setValues(EvalTree& et) const
{
	// set values of constants
	aa.atoms.setValues(et);

	// set values of variables to NaN or to user set values
	double nan = std::numeric_limits<double>::quiet_NaN();
	for (int i = 0; i < aa.atoms.nvar(); i++) {
		const char* ss = aa.atoms.name(i);
		int t = aa.atoms.index(ss);
		if (t >= 0) {
			Tusrvalmap::const_iterator it = user_values.find(t);
			if (it == user_values.end())
				et.set_nulary(t, nan);
			else
				et.set_nulary(t, (*it).second);
		}
	}
}

void AtomAsgnEvaluator::set_user_value(const char* name, double val)
{
	int t = aa.atoms.index(name);
	if (t >= 0) {
		Tusrvalmap::iterator it = user_values.find(t);
		if (it == user_values.end())
			user_values.insert(Tusrvalmap::value_type(t, val));
		else
			(*it).second = val;
	}
}

void AtomAsgnEvaluator::load(int i, double res)
{
	// set the value
	operator[](i) = res;
	// if i-th expression is atom, set its value to this EvalTree
	int t = aa.order[i];
	if (t >= 0)
		etree.set_nulary(t, res);
}

double AtomAsgnEvaluator::get_value(const char* name) const
{
	AtomAssignings::Tvarintmap::const_iterator it = aa.lname2expr.find(name);
	if (it == aa.lname2expr.end())
		return std::numeric_limits<double>::quiet_NaN();
	else
		return operator[]((*it).second);
}


