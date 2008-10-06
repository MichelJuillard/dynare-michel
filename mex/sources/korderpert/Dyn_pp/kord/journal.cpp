#define _W32_FT_OFFSET (116444736000000000LL) 
/*1:*/
#line 6 "./journal.cweb"

#include "journal.h"
#include "kord_exception.h"

#ifndef __MINGW32__
# include <sys/resource.h> 
# include <sys/utsname.h> 
#endif
#include <stdlib.h> 
#include <unistd.h> 
#include <time.h> 

SystemResources _sysres;
#ifdef __MINGW32__
/*16:*/
#line 249 "./journal.cweb"

typedef struct _filetime{
unsigned long dwLowDateTime;
unsigned long dwHighDateTime;
}filetime;

extern"C"{
void __stdcall GetSystemTimeAsFileTime(filetime*);
};

typedef union{
long long ns100;
filetime ft;
}w32_ftv;

void D_gettimeofday(struct timeval*p,struct timezone*tz)
{
w32_ftv _now;
GetSystemTimeAsFileTime(&(_now.ft));
p->tv_usec= (long)((_now.ns100/10LL)%1000000LL);
p->tv_sec= (long)((_now.ns100-_W32_FT_OFFSET)/10000000LL);
return;
}

/*:16*/
#line 20 "./journal.cweb"
;
/*17:*/
#line 282 "./journal.cweb"

#define _SC_PAGESIZE 1
#define _SC_PHYS_PAGES 2
#define _SC_AVPHYS_PAGES 3
#define _SC_NPROCESSORS_ONLN 4

struct Win32MemoryStatus{
unsigned long dwLength;
unsigned long dwMemoryLoad;
unsigned int dwTotalPhys;
unsigned int dwAvailPhys;
unsigned int dwTotalPageFile;
unsigned int dwAvailPageFile;
unsigned int dwTotalVirtual;
unsigned int dwAvailVirtual;
Win32MemoryStatus();
};

extern"C"{
void __stdcall GlobalMemoryStatus(Win32MemoryStatus*);
};

Win32MemoryStatus::Win32MemoryStatus()
{
dwLength= sizeof(Win32MemoryStatus);
GlobalMemoryStatus(this);
}

long sysconf(int name)
{
switch(name){
case _SC_PAGESIZE:
return 1024;
case _SC_PHYS_PAGES:
{
Win32MemoryStatus memstat;
return memstat.dwTotalPhys/1024;
}
case _SC_AVPHYS_PAGES:
{
Win32MemoryStatus memstat;
return memstat.dwAvailPhys/1024;
}
case _SC_NPROCESSORS_ONLN:
return-1;
default:
KORD_RAISE("Not implemented in Win32 sysconf.");
return-1;
}
}

/*:17*/
#line 21 "./journal.cweb"
;
#endif

/*2:*/
#line 40 "./journal.cweb"

SystemResources::SystemResources()
{
D_gettimeofday(&start,NULL);
}


/*:2*/
#line 24 "./journal.cweb"
;
/*3:*/
#line 48 "./journal.cweb"

long int SystemResources::pageSize()
{
return sysconf(_SC_PAGESIZE);
}

/*:3*/
#line 25 "./journal.cweb"
;
/*4:*/
#line 55 "./journal.cweb"

long int SystemResources::physicalPages()
{
return sysconf(_SC_PHYS_PAGES);
}

/*:4*/
#line 26 "./journal.cweb"
;
/*5:*/
#line 62 "./journal.cweb"

long int SystemResources::onlineProcessors()
{
return sysconf(_SC_NPROCESSORS_ONLN);
}

/*:5*/
#line 27 "./journal.cweb"
;
/*6:*/
#line 69 "./journal.cweb"

long int SystemResources::availableMemory()
{
return pageSize()*sysconf(_SC_AVPHYS_PAGES);
}

/*:6*/
#line 28 "./journal.cweb"
;
/*7:*/
#line 78 "./journal.cweb"

void SystemResources::getRUS(double&load_avg,long int&pg_avail,
double&utime,double&stime,double&elapsed,
long int&idrss,long int&majflt)
{
struct timeval now;
D_gettimeofday(&now,NULL);
elapsed= now.tv_sec-start.tv_sec+(now.tv_usec-start.tv_usec)*1.0e-6;

#ifndef __MINGW32__
struct rusage rus;
getrusage(RUSAGE_SELF,&rus);
utime= rus.ru_utime.tv_sec+rus.ru_utime.tv_usec*1.0e-6;
stime= rus.ru_stime.tv_sec+rus.ru_stime.tv_usec*1.0e-6;
idrss= rus.ru_idrss;
majflt= rus.ru_majflt;

getloadavg(&load_avg,1);
#else
utime= -1.0;
stime= -1.0;
idrss= -1;
majflt= -1;
load_avg= -1.0;
#endif
pg_avail= sysconf(_SC_AVPHYS_PAGES);
}

/*:7*/
#line 29 "./journal.cweb"
;
/*8:*/
#line 107 "./journal.cweb"

SystemResourcesFlash::SystemResourcesFlash()
{
_sysres.getRUS(load_avg,pg_avail,utime,stime,
elapsed,idrss,majflt);
}

/*:8*/
#line 30 "./journal.cweb"
;
/*9:*/
#line 115 "./journal.cweb"

void SystemResourcesFlash::diff(const SystemResourcesFlash&pre)
{
utime-= pre.utime;
stime-= pre.stime;
elapsed-= pre.elapsed;
idrss-= pre.idrss;
majflt-= pre.majflt;
}

/*:9*/
#line 31 "./journal.cweb"
;
/*10:*/
#line 126 "./journal.cweb"

JournalRecord&JournalRecord::operator<<(const IntSequence&s)
{
operator<<("[");
for(int i= 0;i<s.size();i++){
operator<<(s[i]);
if(i<s.size()-1)
operator<<(",");
}
operator<<("]");
return*this;
}

/*:10*/
#line 32 "./journal.cweb"
;
/*11:*/
#line 140 "./journal.cweb"

void JournalRecord::writePrefix(const SystemResourcesFlash&f)
{
for(int i= 0;i<MAXLEN;i++)
prefix[i]= ' ';
double mb= 1024*1024;
sprintf(prefix,"%07.6g",f.elapsed);
sprintf(prefix+7,":%c%05d",recChar,ord);
sprintf(prefix+14,":%1.1f",f.load_avg);
sprintf(prefix+18,":%05.4g",f.pg_avail*_sysres.pageSize()/mb);
sprintf(prefix+24,"%s",":      : ");
for(int i= 0;i<2*journal.getDepth();i++)
prefix[i+33]= ' ';
prefix[2*journal.getDepth()+33]= '\0';
}

/*:11*/
#line 33 "./journal.cweb"
;
/*12:*/
#line 157 "./journal.cweb"

void JournalRecordPair::writePrefixForEnd(const SystemResourcesFlash&f)
{
for(int i= 0;i<MAXLEN;i++)
prefix_end[i]= ' ';
double mb= 1024*1024;
SystemResourcesFlash difnow;
difnow.diff(f);
sprintf(prefix_end,"%07.6g",f.elapsed+difnow.elapsed);
sprintf(prefix_end+7,":E%05d",ord);
sprintf(prefix_end+14,":%1.1f",difnow.load_avg);
sprintf(prefix_end+18,":%05.4g",difnow.pg_avail*_sysres.pageSize()/mb);
sprintf(prefix_end+24,":%06.5g",difnow.majflt*_sysres.pageSize()/mb);
sprintf(prefix_end+31,"%s",": ");
for(int i= 0;i<2*journal.getDepth();i++)
prefix_end[i+33]= ' ';
prefix_end[2*journal.getDepth()+33]= '\0';
}

/*:12*/
#line 34 "./journal.cweb"
;
/*13:*/
#line 177 "./journal.cweb"

JournalRecordPair::~JournalRecordPair()
{
journal.decrementDepth();
writePrefixForEnd(flash);
journal<<prefix_end;
journal<<mes;
journal<<endl;
journal.flush();
}

/*:13*/
#line 35 "./journal.cweb"
;
/*14:*/
#line 189 "./journal.cweb"

JournalRecord&endrec(JournalRecord&rec)
{
rec.journal<<rec.prefix;
rec.journal<<rec.mes;
rec.journal<<endl;
rec.journal.flush();
rec.journal.incrementOrd();
return rec;
}

/*:14*/
#line 36 "./journal.cweb"
;
/*15:*/
#line 201 "./journal.cweb"

void Journal::printHeader()
{
(*this)<<"This is Dynare++, Copyright (C) 2004,2005 Michel Juillard, Ondra Kamenik\n";
(*this)<<"Dynare++ comes with ABSOLUTELY NO WARRANTY and is distributed under\n";
(*this)<<"General Public License, see http://www.gnu.org/license/gpl.html\n";
(*this)<<"\n\n";

#ifndef __MINGW32__
utsname info;
uname(&info);
(*this)<<"System info: ";
(*this)<<info.sysname<<" "<<info.release<<" "<<info.version<<" ";
(*this)<<info.machine<<", processors online: "<<_sysres.onlineProcessors();

(*this)<<"\n\nStart time: ";
char ts[100];
time_t curtime= time(NULL);
tm loctime;
localtime_r(&curtime,&loctime);
asctime_r(&loctime,ts);
(*this)<<ts<<"\n";
#else
(*this)<<"System info: (not implemented for MINGW)\n";
(*this)<<"Start time:  (not implemented for MINGW)\n\n";
#endif

(*this)<<"  ------ elapsed time (seconds)                     \n";
(*this)<<"  |       ------ record unique identifier           \n";
(*this)<<"  |       |     ------ load average                 \n";
(*this)<<"  |       |     |    ------ available memory (MB)   \n";
(*this)<<"  |       |     |    |     ------  major faults (MB)\n";
(*this)<<"  |       |     |    |     |                        \n";
(*this)<<"  V       V     V    V     V                        \n";
(*this)<<"\n";
}


/*:15*/
#line 37 "./journal.cweb"
;

/*:1*/
