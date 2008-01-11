/*
 * Copyright (C) 2007-2008 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "linbcg.hh"
//#define DEBUG

LinBCG::LinBCG()
{
 ijat_p=NULL;
}


void LinBCG::Preconditioner(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> IM_i, int* index_vara, int* index_equa, int y_size, double *y, bool print_it, int type, Mat_O_DP &a, Vec_O_INT &indx)
{
  std::map<std::pair<std::pair<int, int>, int> , int>::iterator it4_i;
  std::map<std::pair< int, int >, double>::iterator it4;
  int i, j, k, eq=0, eqp, var, lag, pos_u;
  std::map<std::pair< int, int >, double> IM;
  LinBCG U_linbcg, L_linbcg;
  Vec_INT ijL, ijU;
  Vec_DP sL, sU;
  double d;

  IM.clear();
  it4_i=IM_i.begin();
  int u_size=IM_i.size();
  //mexPrintf("Preconditioner\n");
  double u_pos_u;
  while(it4_i!=IM_i.end())
    {
      eq= it4_i->first.first.first;
      var=it4_i->first.first.second;
      lag=it4_i->first.second;
      pos_u=it4_i->second;
      if(pos_u>=Size)  // to eliminate the residual contain in the first (Size) elements of IM
        {
          //mexPrintf("IM_i[%d, (%d, %d)]=%d val=%f ",eq, var, lag, pos_u, u[pos_u]);
          IM[make_pair(eq,var-lag*Size)]+=u[pos_u/*+(periods-1)*u_size*/];
          //mexPrintf(" IM[%d, %d]=%f\n",eq,var-lag*Size, IM[make_pair(eq,var-lag*Size)]);
        }
      it4_i++;
    }
  int Size_SparseMatrixRow=IM.size();
  //mexPrintf("Size_SparseMatrixRow+1=%d\n",Size_SparseMatrixRow+1);
  int ija_d[Size_SparseMatrixRow+1];
  DP sa_d[Size_SparseMatrixRow+1];
  k=Size;
  ija_d[0]=k+1;
  j=0;
  eqp=-99;
  it4=IM.begin();
  while(it4!=IM.end())
    {
      eq=it4->first.first;
      var=it4->first.second;
      u_pos_u=it4->second;
      //mexPrintf("IM[eq=%d, var=%d]=%f \n",eq, var, u_pos_u);
      if(eq!=eqp)
        ija_d[eq+j*Size]=k+1;
      if(eq==var)
        sa_d[var+j*Size]+=u_pos_u;//u[pos_u+j*u_size];
      else
        {
          k++;
          sa_d[k]+=u_pos_u;//u[pos_u+j*u_size];
          ija_d[k]=var+j*Size;
        }
      eqp=eq;
      it4++;
    }
  eq=Size;//Size*periods;
  ija_d[eq]=k+1;
  sa_d[eq]=0;
  ija_p.Format(ija_d,Size_SparseMatrixRow+1);
  sa_p.Format(sa_d,Size_SparseMatrixRow+1);
  /*sprsprt();
  sprs_sprt();*/
  //mexPrintf("ok0\n");
  a.Format(Size,Size);
  sprsout(a);
  indx.Format(Size);
  //mexPrintf("ok1\n");
  /*Mat_DP aa(0.0,Size*periods,Size*periods);
  for(k=0;k<periods;k++)
    for(i=0;i<Size;i++)
      for(j=0;j<Size;j++)
        aa[i+k*Size][j+k*Size]=a[i][j];
  a=aa;*/
  /*mexPrintf("a before=\n");
  a.Print();*/
  sprsludcmp(a,indx, d);
  /*mexPrintf("a after=\n");
  a.Print();*/
  //mexPrintf("end of preconditioner\n");
  /*sprs_LU(sL, ijL, sU, ijU);

  U_linbcg.sprsin(sU, ijU);

  mexPrintf("U\n");
  U_linbcg.sprs_sprt();
  U_linbcg.sprsprt();

  mexPrintf("ijL[0]=%d\n",ijL[0]);
  L_linbcg.sprsin(sL, ijL);
  mexPrintf("L\n");
  L_linbcg.sprs_sprt();
  L_linbcg.sprsprt();*/
}


void LinBCG::SolveLinear(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> IM_i, int* index_vara, int* index_equa, int y_size, double *y, bool print_it, bool cvg, Mat_DP &a, Vec_IO_INT &indx)
{
  // Solve iteratively the systme A*x=b
  // A is sparse and rewritten using sa_p and ija_p
  mexPrintf("SolveLinear\n");
  mexEvalString("drawnow;");
  std::map<std::pair<std::pair<int, int>, int> , int>::iterator it4_i;
  std::map<std::pair<int, std::pair<int,int> >, int>::iterator it4;
  int i, j, k, eq=0, eqp, var, lag, pos_u;
  std::map<std::pair< int, std::pair<int,int> >, int> IM;
  clock_t t1 = clock();
  IM.clear();
  it4_i=IM_i.begin();
  int u_size=IM_i.size();
  int Size_SparseMatrixRow=0;
  Vec_DP b(Size*periods), x(Size*periods);
  int ITMAX=10*Size*periods;
  double TOL=1e-9;
  int ITOL=1;
  int sub_iter;
  DP err;
  //Vec_DP Err(Size*periods);

  if (iter>0)
    mexPrintf("Sim : %f ms\n",(1000.0*(double(clock())-double(time00)))/double(CLOCKS_PER_SEC));
  if (isnan(res1) || isinf(res1))
    {
      if (slowc_save<1e-8)
        {
          mexPrintf("Dynare cannot improve the simulation\n");
          mexEvalString("drawnow;");
          filename+=" stopped";
          mexEvalString("st=fclose('all');clear all;");
          mexErrMsgTxt(filename.c_str());
        }
      slowc_save/=2;
      mexPrintf("Error: Simulation diverging, trying to correct it using slowc=%f\n",slowc_save);
      for (i=0;i<y_size*(periods+y_kmin);i++)
        y[i]=ya[i]+slowc_save*direction[i];
      iter--;
      return;
    }
  mexPrintf("-----------------------------------\n");
  mexPrintf("      Simulate     iteration° %d     \n",iter+1);
  mexPrintf("      max. error=%.10e       \n",double(max_res));
  mexPrintf("      sqr. error=%.10e       \n",double(res2));
  mexPrintf("      abs. error=%.10e       \n",double(res1));
  mexPrintf("-----------------------------------\n");
  mexPrintf("cvg=%d\n\n",cvg);
  mexEvalString("drawnow;");
  if (cvg)
    return;
  for(j=0;j<periods;j++)
    for(i=0;i<Size;i++)
      {
        //mexPrintf("b[%d]=-u[%d]=%f\n",i+j*Size,i+j*u_size,-u[i+j*u_size]);
        b[i+j*Size]=-u[i+j*u_size];
      }
  for(i=0;i<Size;i++)
    for(j=y_kmin;j<periods+y_kmin;j++)
      {
        //mexPrintf("x[%d]=y[%d]=%f\n",i+(j-y_kmin)*Size, index_vara[i]+j*y_size, y[index_vara[i]+j*y_size]);
        x[i+(j-y_kmin)*Size]=y[index_vara[i]+j*y_size];
      }
  Size_SparseMatrixRow=periods*(u_size-Size);
  while(it4_i!=IM_i.end())
    {
      eq= it4_i->first.first.first;
      var=it4_i->first.first.second;
      lag=it4_i->first.second;
      pos_u=it4_i->second;
      if(pos_u>=Size)  // to eliminate the residual contain in the first (Size) elements of IM
        {
          //mexPrintf("IM[%d, (%d, %d)]=%d\n",eq, var, lag, pos_u);
          IM[make_pair(eq,make_pair(var,lag))]=pos_u;
          for(j=0;j<=min(y_kmax+y_kmin,periods-1);j++)
            {
              //mexPrintf("j=%d\n",j);
              if(j-y_kmin<=0)
                k=j;
              else
                k=periods+j-y_kmin-y_kmax-1;
              if(lag+j<0 || lag+j>min(periods-1,y_kmax+y_kmin))
                {
                  //mexPrintf("eliminate b[eq+(k)*Size=%d]=u[%d]*y[%d]\n",eq+(k)*Size,pos_u+(k)*u_size,var+(k+1)*y_size);
                  b[eq+(k)*Size]-=u[pos_u+k*u_size]*y[index_vara[var-lag*Size]+(k+lag+1)*y_size];
                  Size_SparseMatrixRow--;
                }

            }
        }
      it4_i++;
    }
  //mexPrintf("Size_SparseMatrixRow+1=%d\n",Size_SparseMatrixRow+1);
  int ija_d[Size_SparseMatrixRow+1];
  DP sa_d[Size_SparseMatrixRow+1];
  k=Size*periods;
  ija_d[0]=k+1;
  for(j=0;j<periods;j++)
    {
      eqp=-99;
      it4=IM.begin();
      while(it4!=IM.end())
        {
          eq=it4->first.first;
          var=it4->first.second.first;
          lag=it4->first.second.second;
          pos_u=it4->second;
          /*if(j==1)
            mexPrintf("IM[eq=%d, var=%d, lag=%d]=%d val=%f \n",eq, var, lag, pos_u,u[pos_u+j*u_size]);*/
          if(lag+j>=0 && lag+j<periods)
            {
              if(eq!=eqp)
                {
                  ija_d[eq+j*Size]=k+1;
                  //mexPrintf("ija_d[%d](1)=%d\n",eq+j*Size,k+1);
                }
              if(eq==var)
                {
#ifdef DEBUG
                  sa_d[var+j*Size]=double(pos_u+j*u_size);
#else
                  sa_d[var+j*Size]=u[pos_u+j*u_size];
#endif
                  //mexPrintf("sa_d[%d]=u[%d]=%f\n",var+j*Size,pos_u+j*u_size,u[pos_u+j*u_size]);
                }
              else
                {
                  if(j+lag>=0 && j+lag<periods)
                    {
                      k++;
#ifdef DEBUG
                      sa_d[k]=double(pos_u+j*u_size);
#else
                      sa_d[k]=u[pos_u+j*u_size];
#endif
                      //mexPrintf("sa_d[%d]=u[%d]=%f\n",k,pos_u+j*u_size,u[pos_u+j*u_size]);
                      ija_d[k]=var+j*Size;
                      //mexPrintf("ija_d[%d](3)=%d\n",k,var+j*Size);
                    }
                }
              eqp=eq;
            }
          it4++;
        }
    }
  eq=Size*periods;
  //mexPrintf("eq=%d\n",eq);
  ija_d[eq]=k+1;
  sa_d[eq]=0;
  /*for(j=0;j<Size_SparseMatrixRow+1;j++)
    {
      mexPrintf("%d ",j);
      mexPrintf("ija=%d ",ija_d[j]);
      mexPrintf("sa=%f\n",sa_d[j]);
    }*/
  //ija_p=(Vec_INT*)mxMalloc(sizeof(Vec_INT));
  ija_p.Format(ija_d,Size_SparseMatrixRow+1);
  //ija_p=;
  //sa_p=(Vec_DP*)mxMalloc(sizeof(Vec_DP));
  sa_p.Format(sa_d,Size_SparseMatrixRow+1);
  //sa_p=;
  /*ija_p=new Vec_INT(ija_d,Size_SparseMatrixRow+1);
  sa_p=new Vec_DP(sa_d,Size_SparseMatrixRow+1);*/
  //sprsprt(*sa_p, *ija_p);
  /*atimes(x,Err,0);
  for(i=0;i<Size*periods;i++)
    {
      Err[i]-=b[i];
      //mexPrintf("Err[%d]=%f\n",i,Err[i]);
    }*/
  BiCGStab(b,x,ITOL,TOL,ITMAX,sub_iter,err, a, indx, periods);
  mexPrintf("Estimated error: %f\n",err);
  mexPrintf("Iterations needed: %d\n",sub_iter);
  //mexPrintf("Solution vector:\n");
  for (i=0;i<y_size*(periods+y_kmin);i++)
    ya[i]=y[i];
  slowc_save=slowc;
  //mexPrintf("slowc=%f\n",slowc);
  for(i=0;i<Size;i++)
    for(j=y_kmin;j<periods+y_kmin;j++)
      {
        direction[index_vara[i]+j*y_size]=x[i+(j-y_kmin)*Size]-y[index_vara[i]+j*y_size];
        y[index_vara[i]+j*y_size]+=slowc*direction[index_vara[i]+j*y_size];
        //mexPrintf("x[%d]=y[%d]=%f\n",i+(j-y_kmin)*Size, index_vara[i]+j*y_size, y[index_vara[i]+j*y_size]);
      }
  /*if(ija_p)
    mxFree(ija_p);
  if(sa_p)
    mxFree(sa_p);*/
  time00=clock();
  if (print_it)
    {
      //pctimer_t t2 = pctimer();
      clock_t t2 = clock();
      mexPrintf("(** %f milliseconds **)\n", 1000.0*(double(t2) - double(t1))/double(CLOCKS_PER_SEC));
      mexEvalString("drawnow;");
    }
}




void LinBCG::asolve(Vec_I_DP &b, Vec_O_DP &x, const int itrnsp, Mat_DP &a, Vec_IO_INT &indx, const int periods)
{
	int i;

	int n=b.size();
	//mexPrintf("n=%d a.nrows()=%d a.ncols()=%d\n",n, a.nrows(), a.ncols());
	//DIAG(A)
	//for(i=0;i<n;i++) x[i]=((sa_p)[i] != 0.0 ? b[i]/(sa_p)[i] : b[i]);
	//Identity
	//for(i=0;i<n;i++) x[i]=b[i];
	//LU
	/*mexPrintf("x=\n");
	x.Print();*/

	/*mexPrintf("b=\n");
	Vec_O_DP(b).Print();*/
	/*mexPrintf("x=\n");
	x.Print();*/
	x=b;
	if(itrnsp)
	  {
	    //lubksb_trp(a, indx, x);
	    lubksb_blck_trp(a, indx, x, periods);
	  }
  else
    {
      //lubksb(a, indx, x);
	    lubksb_blck(a, indx, x, periods);
    }
	/*mexPrintf("solution\n");
  x.Print();

  Vec_DP z(n);

  atimes(x,z,itrnsp);

  //for(i=0;i<n;i++)
    //z[i]=z[i]-b[i];

  mexPrintf("check\n");
  z.Print();*/
}

void LinBCG::atimes(Vec_I_DP &x, Vec_O_DP &r, const int itrnsp)
{
	if (itrnsp) sprstx(sa_p,ija_p,x,r);
	else sprsax(sa_p,ija_p,x,r);
}

DP LinBCG::snrm(Vec_I_DP &sx, const int itol)
{
	int i,isamax;
	DP ans;

	int n=sx.size();
	if (itol <= 3) {
		ans = 0.0;
		for (i=0;i<n;i++) ans += sx[i]*sx[i];
		return sqrt(ans);
	} else {
		isamax=0;
		for (i=0;i<n;i++) {
			if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}

void LinBCG::sprsax(Vec_I_DP &sa, Vec_I_INT &ija, Vec_I_DP &x, Vec_O_DP &b)
{
	int i,k;
  //mexPrintf("sprsax\n");
	int n=x.size();
	/*mexPrintf("n=%d\n",n);
	mexPrintf("ija[0]=%d\n",ija[0]);*/
	if (ija[0] != n+1)
		mexPrintf("Error: \n sprsax: mismatched vector and matrix\n");
  /*mexPrintf("ok0a\n");*/
	for (i=0;i<n;i++) {
		b[i]=sa[i]*x[i];
		for (k=ija[i];k<ija[i+1];k++) {
			b[i] += sa[k]*x[ija[k]];
		}
	}
	/*mexPrintf("ok0b\n");*/
}

void LinBCG::sprstx(Vec_I_DP &sa, Vec_I_INT &ija, Vec_I_DP &x, Vec_O_DP &b)
{
	int i,j,k;
  //mexPrintf("sprstx\n");
	int n=x.size();
	if (ija[0] != (n+1))
		mexPrintf("Error: \n mismatched vector and matrix in sprstx\n");
	for (i=0;i<n;i++) b[i]=sa[i]*x[i];
	for (i=0;i<n;i++) {
		for (k=ija[i];k<ija[i+1];k++) {
			j=ija[k];
			b[j] += sa[k]*x[i];
		}
	}
}


void LinBCG::sprsin(Vec_DP s, Vec_INT ij)
{
  sa_p.Format(s);
  ija_p.Format(ij);
}

void LinBCG::sprsin(DP *s, int *ij, int size)
{
  sa_p.Format(size);
  ija_p.Format(size);
  for(int i=0;i<size;i++)
    {
      sa_p[i]=s[i];
      ija_p[i]=ij[i];
    }
}


void LinBCG::sprsin(double *a, int NbRow, int NbCol, const DP thresh)
{
	int i,j,k;
  int n=NbRow;
  for(i=1, j=0;i<NbRow*NbCol;i++)
    {
      //mexPrintf("%d %% %d = %d a[%d]=%f\n",i,NbRow+1, (i)%(NbRow+1),i,a[i]);
      if(fabs(a[i])>=thresh && ((i)%(NbRow+1))!=0)
        j++;
    }
  //mexPrintf("j=%d\n",j);
	ija_p.Format(j+NbRow+1);
	sa_p.Format(j+NbRow+1);
	//mexPrintf("nz=%d\n",j+NbRow+1);
	for(j=0;j<n;j++)
	  sa_p[j]=a[j+j*n];
	ija_p[0]=n+1;
	k=n;
	for(i=0;i<n;i++)
	  {
		  for(j=0;j<n;j++)
		    {
			    if(fabs(a[j*n+i]) >= thresh && i != j)
			      {
				      sa_p[++k]=a[j*n+i];
				      ija_p[k]=j;
			      }
		    }
		  ija_p[i+1]=k+1;
	  }
}

void LinBCG::sprsin(Mat_DP &a, const DP thresh)
{
	int i,j,k=0;

	int n=a.nrows();
	for(i=0;i<n;i++)
	  for(j=0;j<n;j++)
	    if(fabs(a[i][j])>=thresh && j!=i)
	      k++;
	ija_p.Format(k+n+1);
	sa_p.Format(k+n+1);
	for (j=0;j<n;j++) sa_p[j]=a[j][j];
	ija_p[0]=n+1;
	k=n;
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
			if (fabs(a[i][j]) >= thresh && i != j) {
				//if (++k > nmax) nrerror("sprsin: sa and ija too small");
				k++;
				sa_p[k]=a[i][j];
				ija_p[k]=j;
			}
		}
		ija_p[i+1]=k+1;
	}
}


void LinBCG::sprsout(Mat_DP &a)
{
  int i,j,k,l,ke;
  double maxi=0;
  int N=ija_p[0]-1;
  //int total=ija_p[N];
  //mexPrintf("sprsout\n");
  for(i=0;i<N;i++)
    {
      l=ija_p[i];
      k=ija_p[l];
      ke=ija_p[i+1];
      for(j=0;j<N;j++)
        {
          //mexPrintf("i=%d, j=%d\n",i,j);
          if(i==j)
            //mexPrintf(fmt,sa_p[i]);
            a[i][i]=sa_p[i];
          else if(j==k)
            {
              //mexPrintf(fmt,sa_p[l]);
              a[i][j]=sa_p[l];
              if(l+1<ke)
                {
                  l++;
                  k=ija_p[l];
                }
            }
          else
            //mexPrintf(fmt,0.0);
            a[i][j]=0;
        }
      //mexPrintf("\n");
    }
  //mexPrintf("\n");
  //mexPrintf("end of sprsout\n");
}



void LinBCG::sprsprt()
{
  int i,j,k,l,ke;
  double maxi=0;
  int N=ija_p[0]-1;
  int total=ija_p[N];
  char fmt[12];
  memset(fmt,'\0',12);
  string pf;
  for(i=0;i<total;i++)
    if(fabs(sa_p[i])>maxi)
      maxi=fabs(sa_p[i]);
  int size=int(log10(maxi))+1;
  if(size<7)
    {
      itoa(size+8,fmt,10);
      pf= fmt;
      pf.insert(0,"% ");
      pf.append(".6f ");
      pf.copy(fmt,8);
    }
  else
    strcpy(fmt,"%e ");
  for(i=0;i<N;i++)
    {
      l=ija_p[i];
      k=ija_p[l];
      ke=ija_p[i+1];
      for(j=0;j<N;j++)
        {
          if(i==j)
            mexPrintf(fmt,sa_p[i]);
          else if(j==k)
            {
              mexPrintf(fmt,sa_p[l]);
              if(l+1<ke)
                {
                  l++;
                  k=ija_p[l];
                }
            }
          else
            mexPrintf(fmt,0.0);
        }
      mexPrintf("\n");
    }
  mexPrintf("\n");
}


void LinBCG::sprs_sprt()
{
  int i,j,k,l,ke;
  double maxi=0;
  int N=ija_p[0]-1;
  mexPrintf("N=%d\n",N);
  int total=ija_p[N];
  mexPrintf("total=%d\n",total);
  char fmt[12], fmts[12], fmtd[12];
  memset(fmt,'\0',12);
  memset(fmts,'\0',12);
  memset(fmtd,'\0',12);
  string pf, pf1, pf2;
  for(i=0;i<total;i++)
    if(fabs(sa_p[i])>maxi)
      maxi=fabs(sa_p[i]);
  int size=int(log10(maxi))+1;
  if(size<7)
    {
      itoa(size+8,fmt,10);
      pf= fmt;
      pf.insert(0,"% ");
      pf1=pf2=pf;
      pf.append(".6f ");
      pf1.append("s ");
      pf2.append("d ");
      pf.copy(fmt,8);
      pf1.copy(fmts,8);
      pf2.copy(fmtd,8);
    }
  else
    strcpy(fmt,"%e ");
  mexPrintf(fmts,"i");
  mexPrintf(fmts,"ija");
  mexPrintf(fmts,"sa");
  if(ijat_p.size())
    {
      mexPrintf(fmts,"ijat");
      mexPrintf(fmts,"sat");
    }
  mexPrintf("\n");
  for(i=0;i<total;i++)
    {
      mexPrintf(fmtd,i);
      mexPrintf(fmtd,ija_p[i]);
      mexPrintf(fmt,sa_p[i]);
      if(ijat_p.size())
        {
          mexPrintf(fmtd,ijat_p[i]);
          mexPrintf(fmt,sat_p[i]);
        }
      mexPrintf("\n");
    }
  mexPrintf("\n");
}

void LinBCG::BiCG(Vec_I_DP &b, Vec_IO_DP &x, const int itol, const DP tol,
	const int itmax, int &iter, DP &err, Mat_DP &a, Vec_IO_INT &indx, const int periods)
{
	DP ak,akden,bk,bkden=1.0,bknum,bnrm=0,dxnrm,xnrm,zm1nrm,znrm=0;
	const DP EPS=1.0e-14;
	int j;
	int n=b.size();
	//mexPrintf("n=%d\n",n);
	Vec_DP p(n),pp(n),r(n),rr(n),z(n),zz(n);
	iter=0;
	atimes(x,r,0);
	for (j=0;j<n;j++)
	  {
      r[j]=b[j]-r[j];
	    rr[j]=r[j];
	  }
	//atimes(r,rr,0);
	if (itol == 1)
	  {
		  bnrm=snrm(b,itol);
		  asolve(r,z,0, a, indx, periods);
	  }
	else if (itol == 2)
	  {
		  asolve(b,z,0, a, indx, periods);
		  bnrm=snrm(z,itol);
		  asolve(r,z,0, a, indx, periods);
	  }
	else if (itol == 3 || itol == 4)
	  {
		  asolve(b,z,0, a, indx, periods);
		  bnrm=snrm(z,itol);
		  asolve(r,z,0, a, indx, periods);
		  znrm=snrm(z,itol);
	  }
  else
    mexPrintf("Error: \n illegal itol in BiCG\n");
	//cout << fixed << setprecision(6);
	while (iter < itmax)
	  {
		  ++iter;
		  asolve(rr,zz,1, a, indx, periods);
		  for (bknum=0.0,j=0;j<n;j++)
		    bknum += z[j]*rr[j];
		  if (iter == 1)
		    {
			    for (j=0;j<n;j++)
			      {
				      p[j]=z[j];
				      pp[j]=zz[j];
			      }
		    }
		  else
		    {
			    bk=bknum/bkden;
			    for (j=0;j<n;j++)
			      {
				      p[j]=bk*p[j]+z[j];
				      pp[j]=bk*pp[j]+zz[j];
			      }
		    }
		  bkden=bknum;
		  atimes(p,z,0);
		  for(akden=0.0,j=0;j<n;j++)
		    akden += z[j]*pp[j];
		  ak=bknum/akden;
		  atimes(pp,zz,1);
		  for (j=0;j<n;j++)
		    {
			    x[j] += ak*p[j];
			    r[j] -= ak*z[j];
			    rr[j] -= ak*zz[j];
		    }
		  asolve(r,z,0, a, indx, periods);
		  if (itol == 1)
			  err=snrm(r,itol)/bnrm;
		  else if (itol == 2)
			  err=snrm(z,itol)/bnrm;
		  else if (itol == 3 || itol == 4)
		    {
			    zm1nrm=znrm;
			    znrm=snrm(z,itol);
			    if (fabs(zm1nrm-znrm) > EPS*znrm)
			      {
				      dxnrm=fabs(ak)*snrm(p,itol);
				      err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			      }
			    else
			      {
				      err=znrm/bnrm;
				      continue;
			      }
			    xnrm=snrm(x,itol);
			    if (err <= 0.5*xnrm)
			      err /= xnrm;
			    else
			      {
				      err=znrm/bnrm;
				      continue;
			      }
		    }
		  //cout << "iter=" << setw(4) << iter+1 << setw(12) << err << endl;
		  mexPrintf("iter=%d err=%f ak=%f akden=%f\n",iter,err,ak,akden);
		  if (err <= tol)
		    break;
	  }
}

void LinBCG::BiCGStab(Vec_I_DP &b, Vec_IO_DP &x, const int itol, const DP tol,
	const int itmax, int &iter, DP &err, Mat_DP &a, Vec_IO_INT &indx, const int periods)
{
	DP akden,bkden=1.0,bknum,bnrm=0,dxnrm,xnrm,zm1nrm,znrm=0, alpha, omega, beta, omega_num, omega_den, snorm;
	const DP EPS=1.0e-14;
	int j;
	int n=b.size();
	//mexPrintf("n=%d\n",n);
	Vec_DP p(n),v(n),r(n),rr(n),z(n),s(n), t(n), phat(n), shat(n);
	iter=0;
	// Compute the residual
	atimes(x,r,0);
	for (j=0;j<n;j++)
	  {
      r[j]=b[j]-r[j];
	    rr[j]=r[j];
	  }
	//atimes(r,rr,0);
	if (itol == 1)
	  {
		  bnrm=snrm(b,itol);
		  //mexPrintf("asolve(r,z,..)\n");
		  asolve(r,z,0, a, indx, periods);
	  }
	else if (itol == 2)
	  {
		  asolve(b,z,0, a, indx, periods);
		  bnrm=snrm(z,itol);
		  asolve(r,z,0, a, indx, periods);
	  }
	else if (itol == 3 || itol == 4)
	  {
		  asolve(b,z,0, a, indx, periods);
		  bnrm=snrm(z,itol);
		  asolve(r,z,0, a, indx, periods);
		  znrm=snrm(z,itol);
	  }
  else
    mexPrintf("Error: \n illegal itol in BiCGStab\n");
	//cout << fixed << setprecision(6);
	while (iter < itmax)
	  {
		  ++iter;
		  //asolve(rr,zz,1, a, indx);
		  for (bknum=0.0,j=0;j<n;j++)
		    bknum += r[j]*rr[j];
		  if (iter == 1)
		    {
			    for (j=0;j<n;j++)
			      {
			        p[j]=r[j];
			      }
		    }
		  else
		    {
			    beta=bknum/bkden*(alpha/omega);
			    for (j=0;j<n;j++)
			      {
			        p[j]=r[j]+beta*(p[j]-omega*v[j]);
			      }
		    }
      //mexPrintf("asolve(p,phat,..)\n");
      asolve(p,phat,1, a, indx, periods);

      atimes(phat,v,0);
      for(akden=0.0,j=0;j<n;j++)
		    akden += rr[j]*v[j];
      alpha=bknum/akden;
      for(snorm=0.0, j=0;j<n;j++)
        {
          s[j]=r[j]-alpha*v[j];
          snorm+=s[j]*s[j];
        }
      if(sqrt(snorm)<tol)
        {
          for(j=0;j<n;j++)
            x[j]+=alpha*phat[j]+omega*shat[j];
          break;
        }
      //mexPrintf("asolve(s,shat,..)\n");
      asolve(s,shat,1, a, indx, periods);
      atimes(shat,t,0);
      for(omega_num=0.0,omega_den=0.0,j=0;j<n;j++)
        {
          omega_num+=t[j]*s[j];
          omega_den+=t[j]*t[j];
        }
      omega=omega_num/omega_den;
      for(j=0;j<n;j++)
        {
          x[j]+=alpha*phat[j]+omega*shat[j];
          r[j]=s[j]-omega*t[j];
        }
		  bkden=bknum;


		  if (itol == 1)
			  err=snrm(r,itol)/bnrm;
		  else if (itol == 2)
			  err=snrm(z,itol)/bnrm;
		  else if (itol == 3 || itol == 4)
		    {
			    zm1nrm=znrm;
			    znrm=snrm(z,itol);
			    if (fabs(zm1nrm-znrm) > EPS*znrm)
			      {
				      dxnrm=fabs(alpha)*snrm(p,itol);
				      err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			      }
			    else
			      {
				      err=znrm/bnrm;
				      continue;
			      }
			    xnrm=snrm(x,itol);
			    if (err <= 0.5*xnrm)
			      err /= xnrm;
			    else
			      {
				      err=znrm/bnrm;
				      continue;
			      }
		    }
		  //cout << "iter=" << setw(4) << iter+1 << setw(12) << err << endl;
		  mexPrintf("iter=%d err=%f alpha=%f akden=%f\n",iter,err,alpha,akden);
		  mexEvalString("drawnow;");
		  if (err <= tol)
		    break;
	  }
}


void LinBCG::sprs_col_index()
{
  int i, j, k, begin, end;
  int nze=sa_p.size();
  ijat_p.Format(nze);
  sat_p.Format(nze);
  begin=ija_p[0];
  int n=begin-1;
  //mexPrintf("n=%d\n",n);
  map<std::pair<int, int>, double> list;
  map<std::pair<int, int>, double>::iterator list_it;
  for(j=0;j<n;j++)
    sat_p[j]=sa_p[j];
  sat_p[n]=sa_p[n];
  ijat_p[n]=ija_p[n];
  for(i=0;i<n;i++)
    {
      end=ija_p[i+1];
      for(k=begin;k<end;k++)
        {
          //mexPrintf("make_pair(ija_p[k]=%d, i=%d)=%f\n",ija_p[k],i,sa_p[k]);
          list[make_pair(ija_p[k],i)]=sa_p[k];
        }
      begin=end;
    }
  int col=-1;
  int cur_pos=n+1;
  list_it=list.begin();
  while(list_it!=list.end())
    {
      //mexPrintf("list_it->first.first=%d, list_it->first.second=%d, list_it->second=%f\n",list_it->first.first, list_it->first.second, list_it->second);
      while(col!=list_it->first.first)
        {
          //while(col!=list_it->first.first
          col++;
          ijat_p[col]=cur_pos;
        }
      ijat_p[cur_pos]=list_it->first.second;
      sat_p[cur_pos++]=list_it->second;
      list_it++;
    }
}


void LinBCG::sprs_swap_line_copy(map<std::pair<int, int>, double> &list, int &pos_2_blck, int begin, int end)
{
  int i, j, init_pos;
  for(i=begin;i<end;i++)
    {
      init_pos=pos_2_blck;
      list[make_pair(i,init_pos)]=sa_p[i];
      for(j=ija_p[i];j<ija_p[i+1];j++)
        {
          //mexPrintf("list[%d,%d]=%f\n",pos_2_blck,ija_p[j],sa_p[j]);
          list[make_pair(pos_2_blck++,ija_p[j])]=sa_p[j];
        }
    }
}


void LinBCG::sprs_swap_line_exchange(map<std::pair<int, int>, double> &list, int &pos_2_blck, int LS, int LD)
{
  int j;
  bool OK=false, OK1=true;
  int init_pos=pos_2_blck;
  for(j=ija_p[LS];j<ija_p[LS+1];j++)
    {
      if(ija_p[j]==LD)
        {
          OK=true;
          list[make_pair(LD,init_pos)]=sa_p[j];
          //mexPrintf("found list[%d,%d]=%f\n",LD,init_pos,sa_p[j]);
        }
      else
        {
          if(ija_p[j]>=LS && OK1)
            {
              OK1=false;
              if(sa_p[LS]!=0.0)
                {
                  //mexPrintf("list[%d,%d]=%f\n",pos_2_blck,LS,sa_p[LS]);
                  list[make_pair(pos_2_blck++,LS)]=sa_p[LS];
                }
            }
          //mexPrintf("list[%d,%d]=%f\n",pos_2_blck,ija_p[j],sa_p[j]);
          list[make_pair(pos_2_blck++,ija_p[j])]=sa_p[j];
        }
    }
  if(OK1)
    {
      if(sa_p[LS]!=0.0)
        {
          //mexPrintf("list[%d,%d]=%f\n",pos_2_blck,LS,sa_p[LS]);
          list[make_pair(pos_2_blck++,LS)]=sa_p[LS];
        }
    }
  if(!OK)
    {
      list[make_pair(LD,init_pos)]=0.0;
      //mexPrintf("not found list[%d,%d]=%f\n",LD,init_pos,0.0);
    }
}

void LinBCG::sprs_swap_line(int L0, int L1)
{
  int i,j,k, pos_2_blck, init_pos;
  map<std::pair<int, int>, double> list;
  map<std::pair<int, int>, double>::iterator list_it;
  int n=ija_p[0]-1;
  if(L0>L1)
    {
      i=L0;
      L0=L1;
      L1=i;
    }
  else if(L0==L1)
    return;
  else if(L0>=n || L1>=n)
    {
      mexPrintf("out of range in sprs_swap_line (L0=%d and L1=%d)\n",L0, L1);
      mexEvalString("drawnow;");
      mexEvalString("st=fclose('all');clear all;");
      filename+=" stopped";
      mexErrMsgTxt(filename.c_str());
    }
  int nze=sa_p.size();
  bool OK, OK1;
  int nze0=0, nze1=0;
  pos_2_blck=n+1;

  sprs_swap_line_copy(list, pos_2_blck, 0, L0);

  //mexPrintf("L1=>L0\n");

  sprs_swap_line_exchange(list, pos_2_blck, L1, L0);

  sprs_swap_line_copy(list, pos_2_blck, L0+1, L1);

  //mexPrintf("L0=>L1\n");

  sprs_swap_line_exchange(list, pos_2_blck, L0, L1);

  sprs_swap_line_copy(list, pos_2_blck, L1+1, n);

  sa_p.Format(list.size()+1);
  ija_p.Format(list.size()+1);
  //mexPrintf("list.size()+1=%d\n",list.size()+1);
  list_it=list.begin();
  while(list_it!=list.end())
    {
      //mexPrintf("sa_p[list_it->first.first=%d]=list_it->second=%f\n",list_it->first.first,list_it->second);
      sa_p[list_it->first.first]=list_it->second;
      //mexPrintf("ija_p[list_it->first.first=%d]=list_it->first.second=%d\n",list_it->first.first,list_it->first.second);
      ija_p[list_it->first.first]=list_it->first.second;
      list_it++;
    }
  ija_p[n]=list.size()+1;
  sa_p[n]=0;
  sprs_col_index();
}



void LinBCG::sprsludcmp(Mat_DP &a, Vec_O_INT &indx, DP &d)
{
	const DP TINY=1.0e-20;
	int i,imax,j,k, ii, ik, begin_ija, end_ija, begin_ijat, end_ijat, dum1, dum2;
	DP big,dum,sum,temp;
	Vec_DP saa_p(sa_p);
	Vec_INT ijaa_p(ija_p);
  int nze=sa_p.size();
	int n=ija_p[0]-1;
	Vec_DP vv(n);
	d=1.0;

  n=a.nrows();
	//a=(double**)mxMalloc(n*n*sizeof(double));
  //mexPrintf("sprsludcmp\n");

	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
		  {
		    //mexPrintf("%f ",a[i][j]);
			  if ((temp=fabs(a[i][j])) > big) big=temp;
		  }
    //mexPrintf("\n");
		if(big == 0.0)
      {
        mexPrintf("Singular Matrix in routine sprsludcmp\n");
        mexEvalString("st=fclose('all');clear all;");
        mexErrMsgTxt("End of simulate");
      }
		vv[i]=1.0/big;
	}
  /*mexPrintf("First Step\n");
  mexPrintf("a.ncols()=%d\n",a.ncols());
  mexPrintf("a.nrows()=%d\n",a.nrows());*/
	/*for(i=0;i<n;i++)
	  {
	    big=fabs(sa_p[i]);
	    for(j=ija_p[i];j<ija_p[i+1];j++)
	      if(temp=fabs(sa_p[j])>big)
	        big=temp;
      if(big==0.0)
        {
          mexPrintf("Singular Matrix in routine sprsludcmp\n");
          mexErrMsgTxt("End of simulate");
        }
      vv[i]=1.0/big;
	  }*/
  /*sprs_col_index();  // create a col-index compact storage
  for (j=0;j<n;j++)
    {
      //imax=j;
		  for (i=ijat_p[j];i<ijat_p[j+1] && ijat_p[i]<j;i++)
		    {
			    sum=sat_p[i];
          begin_ija=ija_p[i];
          end_ija=ija_p[i+1];
          begin_ijat=ijat_p[i];
          end_ijat=ijat_p[i+1];
          while(begin_ija<end_ija && begin_ijat<end_ijat && begin_ijat<begin_ija)
            {
              if(ija_p[begin_ija]<ijat_p[begin_ijat])
                begin_ija++;
              else if(ija_p[begin_ija]>ijat_p[begin_ijat])
                begin_ijat++;
              else
                sum -= sa_p[begin_ija++]*sat_p[begin_ijat++];
            }
          sat_p[i]=sum;
		    }
      big=0.0;
      for (;i<ijat_p[j+1] ;i++)
		    {
          sum=sat_p[i];
          begin_ija=ija_p[i];
          end_ija=ija_p[i+1];
          begin_ijat=ijat_p[i];
          end_ijat=ijat_p[i+1];
          while(begin_ija<end_ija && ija_p[begin_ija]<j && begin_ijat<end_ijat && ijat_p[begin_ijat]<j)
            {
              if(ija_p[begin_ija]<ijat_p[begin_ijat])
                begin_ija++;
              else if(ija_p[begin_ija]>ijat_p[begin_ijat])
                begin_ijat++;
              else
                sum -= sa_p[begin_ija++]*sat_p[begin_ijat++];
            }
          sat_p[i]=sum;
          if((dum=vv[ija_p[i]]*fabs(sum)) >= big)
            {
				      big=dum;
				      mexPrintf("imax=i=%d\n",i);
				      imax=ija_p[i];
			      }
		    }
     if (j != imax)
       {
         sprs_swap_line(j,imax);
			   //for (k=0;k<n;k++)
			   //  {
			   //    dum=a[imax][k];
				 //    a[imax][k]=a[j][k];
				 //    a[j][k]=dum;
			   //  }
			   d = -d;
			   vv[imax]=vv[j];
       }
     indx[j]=imax;
     if (sa_p[j] == 0.0) sa_p[j]=TINY;
		   if (j != n-1)
  		   {
	  		   dum=1.0/(sat_p[j]);
		  	   for (i=ija_p[j];i<ija_p[j+1];i++)
		  	     if(ija_p[j]>j)
		  	       sa_p[i] *= dum;
		     }
    }
}*/
	for (j=0;j<n;j++)
	  {
		  for (i=0;i<j;i++)
		    {
			    sum=a[i][j];
			    for (k=0;k<i;k++)
			      sum -= a[i][k]*a[k][j];
			    a[i][j]=sum;
		    }
		  big=0.0;
		  //mexPrintf("ok0\n");
		  for (i=j;i<n;i++)
		    {
			    sum=a[i][j];
			    for (k=0;k<j;k++)
			      sum -= a[i][k]*a[k][j];
			    a[i][j]=sum;
			    if ((dum=vv[i]*fabs(sum)) >= big)
			      {
				      big=dum;
				      imax=i;
			      }
		    }
      //mexPrintf("ok1\n");
		  if (j != imax)
		    {
			    for (k=0;k<n;k++)
			      {
				      dum=a[imax][k];
				      a[imax][k]=a[j][k];
				      a[j][k]=dum;
			      }
			    d = -d;
			    vv[imax]=vv[j];
		    }
		  indx[j]=imax;
		  //mexPrintf("ok2\n");
		  if (a[j][j] == 0.0)
		    a[j][j]=TINY;
		  if (j != n-1)
		    {
			    dum=1.0/(a[j][j]);
			    for (i=j+1;i<n;i++)
			      a[i][j] *= dum;
		    }
      //mexPrintf("ok3\n");
	  }
  mexPrintf("end of sprsludcmp\n");
  //sprsin(a,1e-10);
  //mxFree(a);
}



void LinBCG::lubksb_blck(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b, const int periods)
{
	int i,ii,ip,j, k;
	DP sum;
	int n=a.nrows();
	for(k=0;k<periods;k++)
	  {
	    ii=0;
	    for (i=0;i<n;i++)
	      {
    		  ip=indx[i];
	    	  sum=b[ip+k*n];
		      b[ip+k*n]=b[i+k*n];
		      /*if (ii != 0)
			      for (j=ii-1;j<i;j++)
			        sum -= a[i][j]*b[j+k*n];
		      else if (sum != 0.0)
			      ii=i+1;*/
		      for (j=0;j<i;j++)
		        sum -= a[i][j]*b[j+k*n];
		      b[i+k*n]=sum;
	      }
	  }
  for(k=periods-1;k>=0;k--)
    {
	    for (i=n-1;i>=0;i--)
	      {
		      sum=b[i+k*n];
		      for (j=i+1;j<n;j++)
		        sum -= a[i][j]*b[j+k*n];
          if(fabs(a[i][i])<tol)
            {
              mexPrintf("Error: Singularity in lubksb_blck a[%d][%d](max=%d)=%g\n",i,i,n,a[i][i]);
              mexEvalString("st=fclose('all');clear all;");
              mexErrMsgTxt("Exit from Dynare");
            }
		      b[i+k*n]=sum/a[i][i];
	      }
	  }
}

void LinBCG::lubksb_blck_trp(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b, const int periods)
{
	int i,ii,ip,j, k;
	DP sum;
	int n=a.nrows();
	for(k=0;k<periods;k++)
	  {
	    ii=0;
	    for (i=0;i<n;i++)
	      {
    		  ip=indx[i];
	    	  sum=b[ip+k*n];
		      b[ip+k*n]=b[i+k*n];
		      /*if (ii != 0)
			      for (j=ii-1;j<i;j++)
			        sum -= a[j][i]*b[j+k*n];
		      else if (sum != 0.0)
			      ii=i+1;*/
  	      for (j=0;j<i;j++)
			      sum -= a[j][i]*b[j+k*n];
          if(fabs(a[i][i])<tol)
            {
              mexPrintf("Error: Singularity in lubksb_blck_trp a[%d][%d](max=%d)=%g\n",i,i,n,a[i][i]);
              mexEvalString("st=fclose('all');clear all;");
              mexErrMsgTxt("Exit from Dynare");
            }
		      b[i+k*n]=sum/a[i][i];
	      }
	  }
  for(k=periods-1;k>=0;k--)
    {
	    for (i=n-1;i>=0;i--)
	      {
		      sum=b[i+k*n];
		      for (j=i+1;j<n;j++)
		        sum -= a[j][i]*b[j+k*n];
		      b[i+k*n]=sum;
	      }
	  }
}



void LinBCG::lubksb(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b)
{
  // Solving A.x=b
	int i,ii,ip,j;
	DP sum;
	int n=a.nrows();
  ii=0;
  // Solving L.y=b
  // L is a lower triangular matrix with element on the main diagonal equal to 1
  // Implementation of a foreward substitution to solve this system
  for (i=0;i<n;i++)
    {
 		  ip=indx[i];
  	  sum=b[ip];
      b[ip]=b[i];
      if (ii != 0)
	      for (j=ii-1;j<i;j++)
	        sum -= a[i][j]*b[j];
	    else if (sum != 0.0)
			  ii=i+1;
		  b[i]=sum;
    }
  // Solving U.x=y
  // U is an upper triangular matrix with element on the main diagonal is different from 1
  // Implementation of a backward substitution to solve this system
  for (i=n-1;i>=0;i--)
    {
      sum=b[i];
	    for (j=i+1;j<n;j++)
	      sum -= a[i][j]*b[j];
      if(fabs(a[i][i])<tol)
        {
          mexPrintf("Error: Singularity in lubksb a[%d][%d](max=%d)=%g\n",i,i,n,a[i][i]);
          mexEvalString("st=fclose('all');clear all;");
          mexErrMsgTxt("Exit from Dynare");
        }
	    b[i]=sum/a[i][i];
	  }
}


void LinBCG::lubksb_trp(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b)
{
  // Solving A'.x=b <=> (LU)'.x=b <=> U'L'x=b <=> U'.y=b with x solution of L'x=y
	int i,ii,ip,j;
	DP sum;
	int n=a.nrows();
  ii=0;
  // Solving U'.y=b
  // U' is a lower triangular matrix with element on the main diagonal differnet from 1
  // Implementation of a foreward substitution to solve this system
  for (i=0;i<n;i++)
    {
 		  ip=indx[i];
  	  sum=b[ip];
      b[ip]=b[i];
      if (ii != 0)
	      for (j=ii-1;j<i;j++)
	        sum -= a[j][i]*b[j];
	    else if (sum != 0.0)
			  ii=i+1;
      if(fabs(a[i][i])<tol)
        {
          mexPrintf("Error: Singularity in lubksb_trp a[%d][%d](max=%d)=%g\n",i,i,n,a[i][i]);
          mexEvalString("st=fclose('all');clear all;");
          mexErrMsgTxt("Exit from Dynare");
        }
		  b[i]=sum/a[i][i];
    }
  // Solving L'.x=y
  // L' is an upper triangular matrix with element on the main diagonal equal to 1
  // Implementation of a backward substitution to solve this system
  for (i=n-1;i>=0;i--)
    {
      sum=b[i];
	    for (j=i+1;j<n;j++)
	      sum -= a[j][i]*b[j];
	    b[i]=sum;
	  }
}


void LinBCG::sprs_LU(Vec_O_DP &sL, Vec_O_INT &ijL, Vec_O_DP &sU, Vec_O_INT &ijU)
{
  int j,k, nU, nL;
  int n=ija_p[0]-1;
  int L_count, A_count, U_count;
  L_count=A_count=U_count=n+1;
  int nn=ija_p.size();
  Vec_INT indx(n);
  Vec_INT save_ija_p(ija_p);
  Vec_DP save_sa_p(sa_p);
  DP d;

  mexPrintf("sprs_LU\n");

  Mat_DP a(n,n);

  sprsludcmp(a,indx, d);

  Vec_DP res(n), res_save(n);
  for(int i=0;i<n;i++)
    res_save[i]=res[i]=i;

  a.Print();

  //lubksb(a, indx, res, 0);
  lubksb(a, indx, res);

  mexPrintf("solution\n");
  res.Print();

  Vec_DP z(n);

  atimes(res,z,0);

  mexPrintf("check\n");
  z.Print();

  nU=nL=n;
  for(k=0;k<n;k++)
   for(j=ija_p[k];j<ija_p[k+1];j++)
     if(ija_p[j]<k)
        nU++;
     else
        nL++;
  mexPrintf("nU+1=%d, nL+1=%d sa_p.size()=%d\n",nU+1, nL+1, sa_p.size());
  ijU.Format(nU+1);
  sU.Format(nU+1);

  ijL.Format(nL+1);
  //mexPrintf("ijL[0]=%d\n",ijL[0]);
  sL.Format(nL+1);

  for(k=0;k<n;k++)
   {
     //the diagonal elements
     sU[k]=sa_p[k];
     sL[k]=1.0;

     ijU[k]=U_count;
     ijL[k]=L_count;
     for(j=ija_p[k];j<ija_p[k+1];j++)
       {
         if(ija_p[j]<k)
           {
             mexPrintf("U_count=%d, A_count=%d\n",U_count, A_count);
             ijU[U_count]=ija_p[A_count];
             sU[U_count++]=sa_p[A_count++];
           }
         else
           {
             mexPrintf("L_count=%d, A_count=%d\n",L_count, A_count);
             ijL[L_count]=ija_p[A_count];
             sL[L_count++]=sa_p[A_count++];
           }
       }
   }
  sU[n]=0.0;
  sL[n]=0.0;
  ijU[n]=nU+1;
  ijL[n]=nL+1;
  mexPrintf("ijL[0]=%d\n",ijL[0]);
  mexPrintf("end0 of sprs_LU\n");
  sa_p.Format(save_sa_p);
  ija_p.Format(save_ija_p);
  mexPrintf("end of sprs_LU\n");
}

void LinBCG::Initialize(string filename_arg, double res1_arg, double res2_arg, double max_res_arg, double slowc_arg, double *ya_arg, double *direction_arg, int iter_arg)
{
  filename=filename_arg;
  res1=res1_arg;
  res2=res2_arg;
  max_res=max_res_arg;
  slowc=slowc_arg;
  slowc_save=slowc;
  ya=ya_arg;
  direction=direction_arg;
  iter=iter_arg;
}
