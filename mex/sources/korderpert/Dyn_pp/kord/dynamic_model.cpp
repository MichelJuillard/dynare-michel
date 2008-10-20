/*1:*/

#include "dynamic_model.h"

/*2:*/

void NameList::print()const
{
	for(int i= 0;i<getNum();i++)
		printf("%s\n",getName(i));
}

/*:2*/
;
/*3:*/

void NameList::writeMat4(FILE*fd,const char*vname)const
{
	int maxlen= 0;
	for(int i= 0;i<getNum();i++)
		if(maxlen<(int)strlen(getName(i)))
			maxlen= (int)strlen(getName(i));
		
		if(maxlen==0)
			return;
		
		TwoDMatrix m(getNum(),maxlen);
		for(int i= 0;i<getNum();i++)
			for(int j= 0;j<maxlen;j++)
				if(j<(int)strlen(getName(i)))
					m.get(i,j)= (double)(getName(i)[j]);
				else
					m.get(i,j)= (double)(' ');
				
				Mat4Header header(m,vname,"text matrix");
				header.write(fd);
				fwrite(m.getData().base(),sizeof(double),m.nrows()*m.ncols(),fd);
}

/*:3*/
;
/*4:*/

void NameList::writeMat4Indices(FILE*fd,const char*prefix)const
{
	char tmp[100];
	TwoDMatrix aux(1,1);
	for(int i= 0;i<getNum();i++){
		sprintf(tmp,"%s_i_%s",prefix,getName(i));
		aux.get(0,0)= i+1;
		aux.writeMat4(fd,tmp);
	}
}

/*:4*/
;

/*:1*/
