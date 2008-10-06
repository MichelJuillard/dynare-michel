/*1:*/
#line 77 "./sthread.hweb"

#ifndef STHREAD_H
#define STHREAD_H

#ifdef POSIX_THREADS
#include <pthread.h> 
#endif

#include <stdio.h> 
#include <list> 
#include <map> 

namespace sthread{
using namespace std;

class Empty{};
/*2:*/
#line 117 "./sthread.hweb"

template<bool condition,class Then,class Else> 
struct IF{
typedef Then RET;
};

template<class Then,class Else> 
struct IF<false,Then,Else> {
typedef Else RET;
};



/*:2*/
#line 93 "./sthread.hweb"
;
enum{posix,empty};
template<int> class thread_traits;
template<int> class detach_thread;
/*3:*/
#line 134 "./sthread.hweb"

template<int thread_impl> 
class thread{
typedef thread_traits<thread_impl> _Ttraits;
typedef typename _Ttraits::_Tthread _Tthread;
_Tthread th;
public:
virtual~thread(){}
_Tthread&getThreadIden()
{return th;}
const _Tthread&getThreadIden()const
{return th;}
virtual void operator()()= 0;
void run()
{_Ttraits::run(this);}
void detach_run()
{_Ttraits::detach_run(this);}
void exit()
{_Ttraits::exit();}
};

/*:3*/
#line 97 "./sthread.hweb"
;
/*4:*/
#line 160 "./sthread.hweb"

template<int thread_impl> 
class thread_group{
typedef thread_traits<thread_impl> _Ttraits;
typedef thread<thread_impl> _Ctype;
list<_Ctype*> tlist;
typedef typename list<_Ctype*> ::iterator iterator;
public:
static int max_parallel_threads;
void insert(_Ctype*c)
{tlist.push_back(c);}
/*5:*/
#line 180 "./sthread.hweb"

~thread_group()
{
while(!tlist.empty()){
delete tlist.front();
tlist.pop_front();
}
}

/*:5*/
#line 171 "./sthread.hweb"
;
/*7:*/
#line 216 "./sthread.hweb"

void run()
{
int rem= tlist.size();
iterator pfirst= tlist.begin();
while(rem> 2*max_parallel_threads){
pfirst= run_portion(pfirst,max_parallel_threads);
rem-= max_parallel_threads;
}
if(rem> max_parallel_threads){
pfirst= run_portion(pfirst,rem/2);
rem-= rem/2;
}
run_portion(pfirst,rem);
}




/*:7*/
#line 172 "./sthread.hweb"
;
private:
/*6:*/
#line 192 "./sthread.hweb"

iterator run_portion(iterator start,int n)
{
int c= 0;
for(iterator i= start;c<n;++i,c++){
(*i)->run();
}
iterator ret;
c= 0;
for(ret= start;c<n;++ret,c++){
_Ttraits::join(*ret);
}
return ret;
}


/*:6*/
#line 174 "./sthread.hweb"
;
};

/*:4*/
#line 98 "./sthread.hweb"
;
/*8:*/
#line 238 "./sthread.hweb"

template<int thread_impl> 
struct thread_traits{
typedef typename IF<thread_impl==posix,pthread_t,Empty> ::RET _Tthread;
typedef thread<thread_impl> _Ctype;
typedef detach_thread<thread_impl> _Dtype;
static void run(_Ctype*c);
static void detach_run(_Dtype*c);
static void exit();
static void join(_Ctype*c);
};

/*:8*/
#line 99 "./sthread.hweb"
;
/*9:*/
#line 251 "./sthread.hweb"

struct ltmmkey;
typedef pair<const void*,const char*> mmkey;

template<int thread_impl> 
struct mutex_traits{
typedef typename IF<thread_impl==posix,pthread_mutex_t,Empty> ::RET _Tmutex;
typedef map<mmkey,pair<_Tmutex,int> ,ltmmkey> mutex_int_map;
static void init(_Tmutex&m);
static void lock(_Tmutex&m);
static void unlock(_Tmutex&m);
};

/*:9*/
#line 100 "./sthread.hweb"
;
/*10:*/
#line 283 "./sthread.hweb"

struct ltmmkey{
bool operator()(const mmkey&k1,const mmkey&k2)const
{return k1.first<k2.first||
(k1.first==k2.first&&strcmp(k1.second,k2.second)<0);}
};

template<int thread_impl> 
class mutex_map
:public mutex_traits<thread_impl> ::mutex_int_map
{
typedef typename mutex_traits<thread_impl> ::_Tmutex _Tmutex;
typedef mutex_traits<thread_impl> _Mtraits;
typedef pair<_Tmutex,int> mmval;
typedef map<mmkey,mmval,ltmmkey> _Tparent;
typedef typename _Tparent::iterator iterator;
typedef typename _Tparent::value_type _mvtype;
_Tmutex m;
public:
mutex_map()
{_Mtraits::init(m);}
void insert(const void*c,const char*id,const _Tmutex&m)
{_Tparent::insert(_mvtype(mmkey(c,id),mmval(m,0)));}
bool check(const void*c,const char*id)const
{return _Tparent::find(mmkey(c,id))!=_Tparent::end();}
/*11:*/
#line 318 "./sthread.hweb"

mmval*get(const void*c,const char*id)
{
iterator it= _Tparent::find(mmkey(c,id));
if(it==_Tparent::end())
return NULL;
return&((*it).second);
}

/*:11*/
#line 308 "./sthread.hweb"
;
/*12:*/
#line 331 "./sthread.hweb"

void remove(const void*c,const char*id)
{
iterator it= _Tparent::find(mmkey(c,id));
if(it!=_Tparent::end())
erase(it);
}

/*:12*/
#line 309 "./sthread.hweb"
;
void lock_map()
{_Mtraits::lock(m);}
void unlock_map()
{_Mtraits::unlock(m);}

};

/*:10*/
#line 101 "./sthread.hweb"
;
/*13:*/
#line 346 "./sthread.hweb"

template<int thread_impl> 
class synchro{
typedef typename mutex_traits<thread_impl> ::_Tmutex _Tmutex;
typedef mutex_traits<thread_impl> _Mtraits;
public:
typedef mutex_map<thread_impl> mutex_map_t;
private:
const void*caller;
const char*iden;
mutex_map_t&mutmap;
public:
synchro(const void*c,const char*id,mutex_map_t&mmap)
:caller(c),iden(id),mutmap(mmap)
{lock();}
~synchro()
{unlock();}
private:
/*14:*/
#line 373 "./sthread.hweb"

void lock(){
mutmap.lock_map();
if(!mutmap.check(caller,iden)){
_Tmutex mut;
_Mtraits::init(mut);
mutmap.insert(caller,iden,mut);
}
mutmap.get(caller,iden)->second++;
mutmap.unlock_map();
_Mtraits::lock(mutmap.get(caller,iden)->first);
}

/*:14*/
#line 364 "./sthread.hweb"
;
/*15:*/
#line 390 "./sthread.hweb"

void unlock(){
mutmap.lock_map();
if(mutmap.check(caller,iden)){
_Mtraits::unlock(mutmap.get(caller,iden)->first);
mutmap.get(caller,iden)->second--;
if(mutmap.get(caller,iden)->second==0)
mutmap.remove(caller,iden);
}
mutmap.unlock_map();
}

/*:15*/
#line 365 "./sthread.hweb"
;
};

/*:13*/
#line 102 "./sthread.hweb"
;
/*16:*/
#line 405 "./sthread.hweb"

template<int thread_impl> 
struct cond_traits{
typedef typename IF<thread_impl==posix,pthread_cond_t,Empty> ::RET _Tcond;
typedef typename mutex_traits<thread_impl> ::_Tmutex _Tmutex;
static void init(_Tcond&cond);
static void broadcast(_Tcond&cond);
static void wait(_Tcond&cond,_Tmutex&mutex);
static void destroy(_Tcond&cond);
};

/*:16*/
#line 103 "./sthread.hweb"
;
/*17:*/
#line 425 "./sthread.hweb"

template<int thread_impl> 
class condition_counter{
typedef typename mutex_traits<thread_impl> ::_Tmutex _Tmutex;
typedef typename cond_traits<thread_impl> ::_Tcond _Tcond;
int counter;
_Tmutex mut;
_Tcond cond;
bool changed;
public:
/*18:*/
#line 445 "./sthread.hweb"

condition_counter()
:counter(0),changed(true)
{
mutex_traits<thread_impl> ::init(mut);
cond_traits<thread_impl> ::init(cond);
}

/*:18*/
#line 435 "./sthread.hweb"
;
/*19:*/
#line 456 "./sthread.hweb"

~condition_counter()
{
cond_traits<thread_impl> ::destroy(cond);
}

/*:19*/
#line 436 "./sthread.hweb"
;
/*20:*/
#line 465 "./sthread.hweb"

void increase()
{
mutex_traits<thread_impl> ::lock(mut);
counter++;
changed= true;
cond_traits<thread_impl> ::broadcast(cond);
mutex_traits<thread_impl> ::unlock(mut);
}

/*:20*/
#line 437 "./sthread.hweb"
;
/*21:*/
#line 476 "./sthread.hweb"

void decrease()
{
mutex_traits<thread_impl> ::lock(mut);
counter--;
changed= true;
cond_traits<thread_impl> ::broadcast(cond);
mutex_traits<thread_impl> ::unlock(mut);
}

/*:21*/
#line 438 "./sthread.hweb"
;
/*22:*/
#line 490 "./sthread.hweb"

int waitForChange()
{
mutex_traits<thread_impl> ::lock(mut);
if(!changed){
cond_traits<thread_impl> ::wait(cond,mut);
}
changed= false;
int res= counter;
mutex_traits<thread_impl> ::unlock(mut);
return res;
}


/*:22*/
#line 439 "./sthread.hweb"
;
};

/*:17*/
#line 104 "./sthread.hweb"
;
/*23:*/
#line 509 "./sthread.hweb"

template<int thread_impl> 
class detach_thread:public thread<thread_impl> {
public:
condition_counter<thread_impl> *counter;
detach_thread():counter(NULL){}
void installCounter(condition_counter<thread_impl> *c)
{counter= c;}
void run()
{thread_traits<thread_impl> ::detach_run(this);}
};

/*:23*/
#line 105 "./sthread.hweb"
;
/*24:*/
#line 525 "./sthread.hweb"

template<int thread_impl> 
class detach_thread_group{
typedef thread_traits<thread_impl> _Ttraits;
typedef cond_traits<thread_impl> _Ctraits;
typedef detach_thread<thread_impl> _Ctype;
list<_Ctype*> tlist;
typedef typename list<_Ctype*> ::iterator iterator;
condition_counter<thread_impl> counter;
public:
static int max_parallel_threads;
/*25:*/
#line 542 "./sthread.hweb"

void insert(_Ctype*c)
{
tlist.push_back(c);
c->installCounter(&counter);
}

/*:25*/
#line 536 "./sthread.hweb"
;
/*26:*/
#line 550 "./sthread.hweb"

~detach_thread_group()
{
while(!tlist.empty()){
delete tlist.front();
tlist.pop_front();
}
}

/*:26*/
#line 537 "./sthread.hweb"
;
/*27:*/
#line 566 "./sthread.hweb"

void run()
{
int mpt= max_parallel_threads;
iterator it= tlist.begin();
while(it!=tlist.end()){
if(counter.waitForChange()<mpt){
counter.increase();
(*it)->run();
++it;
}
}
while(counter.waitForChange()> 0){}
}


/*:27*/
#line 538 "./sthread.hweb"
;
};

/*:24*/
#line 106 "./sthread.hweb"
;
#ifdef POSIX_THREADS
/*28:*/
#line 586 "./sthread.hweb"

typedef detach_thread<posix> PosixThread;
typedef detach_thread_group<posix> PosixThreadGroup;
typedef synchro<posix> posix_synchro;
class PosixSynchro:public posix_synchro{
public:
PosixSynchro(const void*c,const char*id);
};

#define THREAD sthread::PosixThread
#define THREAD_GROUP sthread::PosixThreadGroup
#define SYNCHRO sthread::PosixSynchro

/*:28*/
#line 108 "./sthread.hweb"
;
#else
/*29:*/
#line 604 "./sthread.hweb"

typedef thread<empty> NoThread;
typedef thread_group<empty> NoThreadGroup;
typedef synchro<empty> no_synchro;
class NoSynchro{
public:
NoSynchro(const void*c,const char*id){}
~NoSynchro(){}
};

#define THREAD sthread::NoThread
#define THREAD_GROUP sthread::NoThreadGroup
#define SYNCHRO sthread::NoSynchro

/*:29*/
#line 110 "./sthread.hweb"
;
#endif
};

#endif

/*:1*/
