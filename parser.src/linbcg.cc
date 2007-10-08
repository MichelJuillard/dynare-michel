#include "linbcg.hh"

LinBCG::LinBCG()
{
}

void LinBCG::Init(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> IM_i, int* index_vara, int* index_equa)
{
  std::map<std::pair<std::pair<int, int> ,int>, int>::iterator it4;
  int Non_Zero_By_Lag[y_kmax+y_kmin+1];
  memset(Non_Zero_By_Lag, 0, sizeof(int)*(y_kmax+y_kmin+1));
  int i, j, k, nz_diag_count, eq, eqp, eq_count, var, lag, pos_u;
  int *inv_index_equa;
  inv_index_equa=(int*)mxMalloc(Size*sizeof(int));
  std::map<std::pair<std::pair<int, int> ,int>, int> IM;
  IM.clear();
  it4=IM_i.begin();
  while(it4!=IM_i.end())
    {
      eq= it4->first.first.first;
      var=it4->first.first.second;
      lag=it4->first.second;
      pos_u=it4->second;
      if(var<Size*periods+y_kmax+y_kmin)
        IM[make_pair(make_pair(eq, var), lag)]= pos_u;
      it4++;
    }
  
  for(j=0;j<Size;j++)
    inv_index_equa[index_equa[j]]=j;
  it4=IM.begin();
  // count the non zero diagonal elements
  j=0;
  nz_diag_count=0;
  eqp=it4->first.first.first;
  mexPrintf("IM.size()=%d\n",IM.size());
  for(i=0;i<Size;i++)
    {
      mexPrintf("index_vara[%d]=%d index_equa[%d]=%d\n",i,index_vara[i],i,index_equa[i]);
    }
  
  while (it4!=IM.end())
    {
      eq= it4->first.first.first;
      // Assumption no Row is completly null
      //if(eq!=eqp)
      j=inv_index_equa[index_equa[eq]];
      var=it4->first.first.second;
      lag=it4->first.second;
      mexPrintf("eq=%d var=%d lag=%d \n",eq, var, lag);
      if(!lag)
        {
          if(var==index_vara[j] && eq==index_equa[j])
            nz_diag_count++;
          else
            Non_Zero_By_Lag[y_kmin+lag]++;
        }
      else
       Non_Zero_By_Lag[y_kmin+lag]++;
      it4++;
    }
  int Size_SparseMatrixRow=Size*periods;
  mexPrintf("En lag 0 Size*periods=%d\n",Size*periods);
  for(lag=-y_kmin;lag<0;lag++)
    {
      Size_SparseMatrixRow+=max(periods+lag,0)*Non_Zero_By_Lag[y_kmin+lag];
      mexPrintf("En lag %d max(periods+lag,0)*Non_Zero_By_Lag[y_kmin+lag]=%d\n",lag,max(periods+lag,0)*Non_Zero_By_Lag[y_kmin+lag]);
    }
  Size_SparseMatrixRow+=periods*Non_Zero_By_Lag[y_kmin];
  for(lag=1;lag<=y_kmax;lag++)
    {
      Size_SparseMatrixRow+=max(periods-lag,0)*Non_Zero_By_Lag[y_kmin+lag];
      mexPrintf("En lag %d max(periods-lag,0)*Non_Zero_By_Lag[y_kmin+lag]=%d\n",lag,max(periods-lag,0)*Non_Zero_By_Lag[y_kmin+lag]);
    }
  mexPrintf("Size_SparseMatrixRow=%d periods=%d\n",Size_SparseMatrixRow,periods);
  //ija_p=new Vec_INT(Size_SparseMatrixRow);
  //sa_p=new Vec_DP(Size_SparseMatrixRow);
  Vec_INT ija_p(Size_SparseMatrixRow);
  Vec_DP sa_p(Size_SparseMatrixRow);
  int *cumulate_past_block;
  cumulate_past_block=(int*)mxMalloc((periods+1)*sizeof(int));
  memset(cumulate_past_block, 0, sizeof(int)*(periods+1)*Size);  
  mexPrintf("OK-1\n");
  cumulate_past_block[0]=0;
  for(j=1;j<=periods;j++)
    {
      for(i=-y_kmin;i<=y_kmax;i++)
        {
          if(j+i>0 && j+i<periods)
            {
              mexPrintf("cumulate_past_block[%d]+=Non_Zero_By_Lag[%d]\n",j,i+y_kmin);
              cumulate_past_block[j]+=Non_Zero_By_Lag[i+y_kmin];
            }
        }
    }
  for(j=0;j<=periods;j++)
    mexPrintf("cumulate_past_block[%d]=%d\n",j,cumulate_past_block[j]);
  int IM_Size=IM.size();
  k=Size*periods;
  mexPrintf("OK00\n");
  ija_p[1]=Size*periods+2;
  //ija_p[Size*periods+1]=Size_SparseMatrixRow+1;
  eq_count=0;eqp=index_equa[0];
  mexPrintf("OK0\n");
  it4=IM.begin();
  while (it4!=IM.end())
    {
      eq= it4->first.first.first;
      // Assumption no Row is completly null
      /*if(eq!=eqp)
        j++;*/
      var=it4->first.first.second;
      lag=it4->first.second;
      pos_u=it4->second;
      //mexPrintf("eq=%d var=%d lag=%d \n",eq, var, lag);
      if(!lag)
        {
          mexPrintf("var=%d index_vara[%d]=%d eq=%d index_equa[%d]=%d\n",var,eq_count,index_vara[eq_count],eq, eq_count, index_equa[eq_count]);
          if(var==index_vara[eq_count] && eq==index_equa[eq_count])
            //nz_diag_count++;
            {
              //Diagonal elements
              for(j=0;j<periods;j++)
                {
                  sa_p[var+j*Size]=u[pos_u+(j+y_kmin)*IM_Size];
                  mexPrintf("sa_p(d)[%d]=u[%d]=%f\n",var+j*Size,pos_u+(j+y_kmin)*IM_Size,sa_p[var+j*Size]);
                }
            }
          else
            {
              k++;
              for(j=0;j<periods;j++)
                {
                  sa_p[cumulate_past_block[j]+k]=u[pos_u+(j+y_kmin)*IM_Size];
                  mexPrintf("sa_p(hd)[%d]=u[%d]=%f\n",cumulate_past_block[j]+k,pos_u+(j+y_kmin)*IM_Size,sa_p[cumulate_past_block[j]+k]);
                  ija_p[cumulate_past_block[j]+k+1]=Size*j+var+1;
                  mexPrintf("ija_p[%d]=%d\n",cumulate_past_block[j]+k+1,ija_p[cumulate_past_block[j]+k+1]);
                }
            }
        }
      else
        {
          k++;
          for(j=0;j<periods;j++)
            if(j+lag>0 && j+lag<periods)
              {
                //sa_p[cumulate_past_block[j]+k]=u[pos_u+j*IM_Size];
                //ija_p[cumulate_past_block[j]+k]=var+1;
                sa_p[cumulate_past_block[j]+k]=u[pos_u+(j+y_kmin)*IM_Size];
                mexPrintf("sa_p(hdl)[%d]=u[%d]=%f\n",cumulate_past_block[j]+k,pos_u+(j+y_kmin)*IM_Size,sa_p[cumulate_past_block[j]+k]);
                ija_p[cumulate_past_block[j]+k+1]=Size*j+var+1;
                mexPrintf("ija_p[%d]=%d\n",cumulate_past_block[j]+k+1,ija_p[cumulate_past_block[j]+k+1]);
              }
        }
      it4++;
      if(eq!=eqp)
        {
          eq_count=inv_index_equa[eq];
          for(j=0;j<periods;j++)
            {
              ija_p[j*Size+eq_count+1]=cumulate_past_block[j]+k;
              mexPrintf("ija_p[%d]=%d\n",j*Size+eq_count+1,ija_p[j*Size+eq_count]);
            }
        }
      eqp=eq;
    }
  mexPrintf("OK1\n");
  for(j=1;j<=Size_SparseMatrixRow;j++)
    {
      mexPrintf("%d ",j);
      mexPrintf("ija=%d ",ija_p[j]);
      mexPrintf("sa=%f\n",sa_p[j]);
    }
  mexErrMsgTxt("Exit from Dynare");
}


void LinBCG::asolve(Vec_I_DP &b, Vec_O_DP &x, const int itrnsp)
{
	int i;

	int n=b.size();
	for(i=0;i<n;i++) x[i]=((*sa_p)[i] != 0.0 ? b[i]/(*sa_p)[i] : b[i]);
}

void LinBCG::atimes(Vec_I_DP &x, Vec_O_DP &r, const int itrnsp)
{
	if (itrnsp) sprstx(*sa_p,*ija_p,x,r);
	else sprsax(*sa_p,*ija_p,x,r);
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

	int n=x.size();
	if (ija[0] != n+1)
		mexPrintf("Error: \n sprsax: mismatched vector and matrix\n");
	for (i=0;i<n;i++) {
		b[i]=sa[i]*x[i];
		for (k=ija[i];k<ija[i+1];k++) {
			b[i] += sa[k]*x[ija[k]];
		}
	}
}

void LinBCG::sprstx(Vec_I_DP &sa, Vec_I_INT &ija, Vec_I_DP &x, Vec_O_DP &b)
{
	int i,j,k;

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

void LinBCG::linbcg(Vec_I_DP &b, Vec_IO_DP &x, const int itol, const DP tol,
	const int itmax, int &iter, DP &err)
{
	DP ak,akden,bk,bkden=1.0,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	const DP EPS=1.0e-14;
	int j;

	int n=b.size();
	Vec_DP p(n),pp(n),r(n),rr(n),z(n),zz(n);
	iter=0;
	atimes(x,r,0);
	for (j=0;j<n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	//atimes(r,rr,0);
	if (itol == 1) {
		bnrm=snrm(b,itol);
		asolve(r,z,0);
	}
	else if (itol == 2) {
		asolve(b,z,0);
		bnrm=snrm(z,itol);
		asolve(r,z,0);
	}
	else if (itol == 3 || itol == 4) {
		asolve(b,z,0);
		bnrm=snrm(z,itol);
		asolve(r,z,0);
		znrm=snrm(z,itol);
	} else mexPrintf("Error: \n illegal itol in linbcg\n");
	cout << fixed << setprecision(6);
	while (iter < itmax) {
		++iter;
		asolve(rr,zz,1);
		for (bknum=0.0,j=0;j<n;j++) bknum += z[j]*rr[j];
		if (iter == 1) {
			for (j=0;j<n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		} else {
			bk=bknum/bkden;
			for (j=0;j<n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes(p,z,0);
		for (akden=0.0,j=0;j<n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		atimes(pp,zz,1);
		for (j=0;j<n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve(r,z,0);
		if (itol == 1)
			err=snrm(r,itol)/bnrm;
		else if (itol == 2)
			err=snrm(z,itol)/bnrm;
		else if (itol == 3 || itol == 4) {
			zm1nrm=znrm;
			znrm=snrm(z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*snrm(p,itol);
				err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				err=znrm/bnrm;
				continue;
			}
			xnrm=snrm(x,itol);
			if (err <= 0.5*xnrm) err /= xnrm;
			else {
				err=znrm/bnrm;
				continue;
			}
		}
		cout << "iter=" << setw(4) << iter+1 << setw(12) << err << endl;
		if (err <= tol) break;
	}
}

