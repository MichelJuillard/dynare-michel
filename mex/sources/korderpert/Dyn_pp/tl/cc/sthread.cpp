/*1:*/
#line 9 "./sthread.cweb"

#include "sthread.h"

#ifdef POSIX_THREADS
namespace sthread{
template<> 
int thread_group<posix> ::max_parallel_threads= 2;
template<> 
int detach_thread_group<posix> ::max_parallel_threads= 2;
/*2:*/
#line 31 "./sthread.cweb"

/*3:*/
#line 40 "./sthread.cweb"

void*posix_thread_function(void*c);
template<> 
void thread_traits<posix> ::run(_Ctype*c)
{
pthread_create(&(c->getThreadIden()),NULL,posix_thread_function,(void*)c);
}

void*posix_detach_thread_function(void*c);

template<> 
void thread_traits<posix> ::detach_run(_Dtype*c)
{
pthread_attr_t attr;
pthread_attr_init(&attr);
pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
pthread_create(&(c->getThreadIden()),&attr,posix_detach_thread_function,(void*)c);
pthread_attr_destroy(&attr);
}


template<> 
void thread_traits<posix> ::exit()
{
pthread_exit(NULL);
}

template<> 
void thread_traits<posix> ::join(_Ctype*c)
{
pthread_join(c->getThreadIden(),NULL);
}

/*:3*/
#line 32 "./sthread.cweb"
;
/*4:*/
#line 74 "./sthread.cweb"

template<> 
void mutex_traits<posix> ::init(pthread_mutex_t&m)
{
pthread_mutex_init(&m,NULL);
}

template<> 
void mutex_traits<posix> ::lock(pthread_mutex_t&m)
{
pthread_mutex_lock(&m);
}

template<> 
void mutex_traits<posix> ::unlock(pthread_mutex_t&m)
{
pthread_mutex_unlock(&m);
}

/*:4*/
#line 33 "./sthread.cweb"
;
/*5:*/
#line 94 "./sthread.cweb"

template<> 
void cond_traits<posix> ::init(_Tcond&cond)
{
pthread_cond_init(&cond,NULL);
}

template<> 
void cond_traits<posix> ::broadcast(_Tcond&cond)
{
pthread_cond_broadcast(&cond);
}

template<> 
void cond_traits<posix> ::wait(_Tcond&cond,_Tmutex&mutex)
{
pthread_cond_wait(&cond,&mutex);
}

template<> 
void cond_traits<posix> ::destroy(_Tcond&cond)
{
pthread_cond_destroy(&cond);
}


/*:5*/
#line 34 "./sthread.cweb"
;
/*6:*/
#line 123 "./sthread.cweb"

static posix_synchro::mutex_map_t posix_mm;

PosixSynchro::PosixSynchro(const void*c,const char*id)
:posix_synchro(c,id,posix_mm){}

/*:6*/
#line 35 "./sthread.cweb"
;
/*7:*/
#line 131 "./sthread.cweb"

void*posix_thread_function(void*c)
{
thread_traits<posix> ::_Ctype*ct= 
(thread_traits<posix> ::_Ctype*)c;
try{
ct->operator()();
}catch(...){
ct->exit();
}
return NULL;
}

/*:7*/
#line 36 "./sthread.cweb"
;
/*8:*/
#line 145 "./sthread.cweb"

void*posix_detach_thread_function(void*c)
{
thread_traits<posix> ::_Dtype*ct= 
(thread_traits<posix> ::_Dtype*)c;
condition_counter<posix> *counter= ct->counter;
try{
ct->operator()();
}catch(...){
ct->exit();
}
if(counter)
counter->decrease();
return NULL;
}


/*:8*/
#line 37 "./sthread.cweb"
;

/*:2*/
#line 18 "./sthread.cweb"
;
}
#else
namespace sthread{
template<> 
int thread_group<empty> ::max_parallel_threads= 1;
template<> 
int detach_thread_group<empty> ::max_parallel_threads= 1;
/*9:*/
#line 166 "./sthread.cweb"

template<> 
void thread_traits<empty> ::run(_Ctype*c)
{
c->operator()();
}
template<> 
void thread_traits<empty> ::detach_run(_Dtype*c)
{
c->operator()();
}

template<> 
void thread_traits<empty> ::exit()
{
}

template<> 
void thread_traits<empty> ::join(_Ctype*c)
{
}

template<> 
void mutex_traits<empty> ::init(Empty&m)
{
}

template<> 
void mutex_traits<empty> ::lock(Empty&m)
{
}

template<> 
void mutex_traits<empty> ::unlock(Empty&m)
{
}

template<> 
void cond_traits<empty> ::init(_Tcond&cond)
{
}

template<> 
void cond_traits<empty> ::broadcast(_Tcond&cond)
{
}

template<> 
void cond_traits<empty> ::wait(_Tcond&cond,_Tmutex&mutex)
{
}

template<> 
void cond_traits<empty> ::destroy(_Tcond&cond)
{
}

/*:9*/
#line 26 "./sthread.cweb"
;
}
#endif

/*:1*/
