/*1:*/

#include "equivalence.h"
#include "permutation.h"
#include "tl_exception.h"

#include <string.h> 

/*2:*/

/*6:*/

int OrdSequence::operator[](int i)const
{
	TL_RAISE_IF((i<0||i>=length()),
		"Index out of range in OrdSequence::operator[]");
	return data[i];
}

/*:6*/
;
/*7:*/

bool OrdSequence::operator<(const OrdSequence&s)const
{
	double ta= average();
	double sa= s.average();
	return(ta<sa||((ta==sa)&&(operator[](0)> s[0])));
}

/*:7*/
;
/*8:*/

bool OrdSequence::operator==(const OrdSequence&s)const
{
	if(length()!=s.length())
		return false;
	
	int i= 0;
	while(i<length()&&operator[](i)==s[i])
		i++;
	
	return(i==length());
}


/*:8*/
;
/*9:*/

void OrdSequence::add(int i)
{
	vector<int> ::iterator vit= data.begin();
	while(vit!=data.end()&&*vit<i)
		++vit;
	if(vit!=data.end()&&*vit==i)
		return;
	data.insert(vit,i);
}

void OrdSequence::add(const OrdSequence&s)
{
	vector<int> ::const_iterator vit= s.data.begin();
	while(vit!=s.data.end()){
		add(*vit);
		++vit;
	}
}

/*:9*/
;
/*10:*/

bool OrdSequence::has(int i)const
{
	vector<int> ::const_iterator vit= data.begin();
	while(vit!=data.end()){
		if(*vit==i)
			return true;
		++vit;
	}
	return false;
}

/*:10*/
;
/*11:*/

double OrdSequence::average()const
{
	double res= 0;
	for(unsigned int i= 0;i<data.size();i++)
		res+= data[i];
	TL_RAISE_IF(data.size()==0,
		"Attempt to take average of empty class in OrdSequence::average");
	return res/data.size();
}

/*:11*/
;
/*12:*/

void OrdSequence::print(const char*prefix)const
{
	printf("%s",prefix);
	for(unsigned int i= 0;i<data.size();i++)
		printf("%d ",data[i]);
	printf("\n");
}

/*:12*/
;

/*:2*/
;
/*3:*/

/*13:*/

Equivalence::Equivalence(int num)
:n(num)
{
	for(int i= 0;i<num;i++){
		OrdSequence s;
		s.add(i);
		classes.push_back(s);
	}
}

Equivalence::Equivalence(int num,const char*dummy)
:n(num)
{
	OrdSequence s;
	for(int i= 0;i<num;i++)
		s.add(i);
	classes.push_back(s);
}

/*:13*/
;
/*14:*/

Equivalence::Equivalence(const Equivalence&e)
:n(e.n),
classes(e.classes)
{
}

Equivalence::Equivalence(const Equivalence&e,int i1,int i2)
:n(e.n),
classes(e.classes)
{
	seqit s1= find(i1);
	seqit s2= find(i2);
	if(s1!=s2){
		OrdSequence ns(*s1);
		ns.add(*s2);
		classes.erase(s1);
		classes.erase(s2);
		insert(ns);
	}
}

/*:14*/
;
/*17:*/

Equivalence::const_seqit Equivalence::findHaving(int i)const
{
	const_seqit si= classes.begin();
	while(si!=classes.end()){
		if((*si).has(i))
			return si;
		++si;
	}
	TL_RAISE_IF(si==classes.end(),
		"Couldn't find equivalence class in Equivalence::findHaving");
	return si;
}

Equivalence::seqit Equivalence::findHaving(int i)
{
	seqit si= classes.begin();
	while(si!=classes.end()){
		if((*si).has(i))
			return si;
		++si;
	}
	TL_RAISE_IF(si==classes.end(),
		"Couldn't find equivalence class in Equivalence::findHaving");
	return si;
}


/*:17*/
;
/*18:*/

Equivalence::const_seqit Equivalence::find(int j)const
{
	const_seqit si= classes.begin();
	int i= 0;
	while(si!=classes.end()&&i<j){
		++si;
		i++;
	}
	TL_RAISE_IF(si==classes.end(),
		"Couldn't find equivalence class in Equivalence::find");
	return si;
}

Equivalence::seqit Equivalence::find(int j)
{
	seqit si= classes.begin();
	int i= 0;
	while(si!=classes.end()&&i<j){
		++si;
		i++;
	}
	TL_RAISE_IF(si==classes.end(),
		"Couldn't find equivalence class in Equivalence::find");
	return si;
}


/*:18*/
;
/*19:*/

void Equivalence::insert(const OrdSequence&s)
{
	seqit si= classes.begin();
	while(si!=classes.end()&&*si<s)
		++si;
	classes.insert(si,s);
}

/*:19*/
;
/*15:*/

const Equivalence&Equivalence::operator= (const Equivalence&e)
{
	classes.clear();
	n= e.n;
	classes= e.classes;
	return*this;
}

/*:15*/
;
/*16:*/

bool Equivalence::operator==(const Equivalence&e)const
{
	if(!std::operator==(classes,e.classes))
		return false;
	
	if(n!=e.n)
		return false;
	
	return true;
}


/*:16*/
;
/*20:*/

void Equivalence::trace(IntSequence&out,int num)const
{
	int i= 0;
	int nc= 0;
	for(const_seqit it= begin();it!=end()&&nc<num;++it,++nc)
		for(int j= 0;j<(*it).length();j++,i++){
			TL_RAISE_IF(i>=out.size(),
				"Wrong size of output sequence in Equivalence::trace");
			out[i]= (*it)[j];
		}
}

/*:20*/
;
/*21:*/

void Equivalence::trace(IntSequence&out,const Permutation&per)const
{
	TL_RAISE_IF(out.size()!=n,
		"Wrong size of output sequence in Equivalence::trace");
	TL_RAISE_IF(per.size()!=numClasses(),
		"Wrong permutation for permuted Equivalence::trace");
	int i= 0;
	for(int iclass= 0;iclass<numClasses();iclass++){
		const_seqit itper= find(per.getMap()[iclass]);
		for(int j= 0;j<(*itper).length();j++,i++)
			out[i]= (*itper)[j];
	}
}


/*:21*/
;
/*22:*/

void Equivalence::print(const char*prefix)const
{
	int i= 0;
	for(const_seqit it= classes.begin();
	it!=classes.end();
	++it,i++){
		printf("%sclass %d: ",prefix,i);
		(*it).print("");
	}
}

/*:22*/
;

/*:3*/
;
/*4:*/

/*23:*/

EquivalenceSet::EquivalenceSet(int num)
:n(num),
equis()
{
	list<Equivalence> added;
	Equivalence first(n);
	equis.push_back(first);
	addParents(first,added);
	while(!added.empty()){
		addParents(added.front(),added);
		added.pop_front();
	}
	if(n> 1){
		Equivalence last(n,"");
		equis.push_back(last);
	}
}

/*:23*/
;
/*24:*/

bool EquivalenceSet::has(const Equivalence&e)const
{
	list<Equivalence> ::const_reverse_iterator rit= equis.rbegin();
	while(rit!=equis.rend()&&*rit!=e)
		++rit;
	if(rit!=equis.rend())
		return true;
	return false;
}

/*:24*/
;
/*25:*/

void EquivalenceSet::addParents(const Equivalence&e,
								list<Equivalence> &added)
{
	if(e.numClasses()==2||e.numClasses()==1)
		return;
	
	for(int i1= 0;i1<e.numClasses();i1++)
		for(int i2= i1+1;i2<e.numClasses();i2++){
			Equivalence ns(e,i1,i2);
			if(!has(ns)){
				added.push_back(ns);
				equis.push_back(ns);
			}
		}
}

/*:25*/
;
/*26:*/

void EquivalenceSet::print(const char*prefix)const
{
	char tmp[100];
	strcpy(tmp,prefix);
	strcat(tmp,"    ");
	int i= 0;
	for(list<Equivalence> ::const_iterator it= equis.begin();
	it!=equis.end();
	++it,i++){
		printf("%sequivalence %d:(classes %d)\n",prefix,i,(*it).numClasses());
		(*it).print(tmp);
	}
}

/*:26*/
;

/*:4*/
;
/*5:*/

/*27:*/

EquivalenceBundle::EquivalenceBundle(int nmax)
{
	nmax= max(nmax,1);
	generateUpTo(nmax);
}

/*:27*/
;
/*28:*/

EquivalenceBundle::~EquivalenceBundle()
{
	for(unsigned int i= 0;i<bundle.size();i++)
		delete bundle[i];
}

/*:28*/
;
/*29:*/

const EquivalenceSet&EquivalenceBundle::get(int n)const
{
	if(n> (int)(bundle.size())||n<1){
		TL_RAISE("Equivalence set not found in EquivalenceBundle::get");
		return*(bundle[0]);
	}else{
		return*(bundle[n-1]);
	}
}

/*:29*/
;
/*30:*/

void EquivalenceBundle::generateUpTo(int nmax)
{
	int curmax= bundle.size();
	for(int i= curmax+1;i<=nmax;i++)
		bundle.push_back(new EquivalenceSet(i));
}


/*:30*/
;


/*:5*/
;

/*:1*/
