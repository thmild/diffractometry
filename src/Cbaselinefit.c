#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Print.h>

# define TOL 1e-16

void extremfnpreg(double *fn,
		  int *pkloc,
		  int *npks,
		  int *knts,
		  int *exvalindl,
		  int *exvalindr,
		  int *minmax)
{
  int i,j,mima;
  j=0;
  if(fn[pkloc[0]]>fn[0])
    mima=1;
  else mima=-1;
  for(i=0;i<*npks;i++){
    do{
      j++;
    }while(knts[j]<pkloc[i]);
    exvalindl[i]=knts[j-1];
    exvalindr[i]=knts[j];
    minmax[i]=mima;
    if(mima==1)
      mima=-1;
    else
      mima=1;
  }
}


int sortfunction (const void *a, const void *b)
{

  if (*(double *)a<*(double *)b)
    return 1;
  else
    return -1;
}



void exber_maxwdth(double * wss,   /* Ergebnis des weighted smoothing splines */
		   double *diffwss, /* Differenzen von wss */
		   double *gr,      /* Grenze: abs(mad(diffwss)) */
		   int *n,         /* Länge der Daten */
		   int *minmax,    /* minmax aus extrempmreg */
		   int *nex,       /* Anzahl der Extremwerte */
		   int *indl,      /* Vektor exvalindl aus extrempmreg mit null davor und dahinter */
		   int *indr,      /* Vektor exvalindr aus extrempmreg mit null davor und dahinter */
		   int *exent,      /* Ergebnis */
		   double *maxwdth, /* maximale Peakbreite */
		   double *x,      /* x-Koordinaten */
		   double *xex     /* x-Koordinaten der Extremwerte */
	   )
{
  int i,j, kn, k;
  double *diffsort;
  double halfmax;
  halfmax = *maxwdth/2.0;
  for(j=1;j<(*nex+1);j++){
    if(minmax[j]==1){
      i=-1;
      do{
	i+=1;
	if((indl[j]-i)==(indr[j-1]-1)) break;
	if((xex[j-1]-x[indl[j]-i]) > halfmax) break;
      }while(diffwss[indl[j]-i-1]<0.0);
      indl[j]-=i;
      i=-1;
      do{
	i+=1;
	if((indr[j]+i)==(indl[j+1]+1)) break;
	if(x[indr[j]+i]-xex[j-1]>halfmax) break;
      }while(diffwss[indr[j]+i-1]>0.0);
      indr[j]+=i;
    }
    else{
      i=-1;
      do{
	i+=1;
	if((indl[j]-i)==(indr[j-1]-1)) break;
      }while(diffwss[indl[j]-i-1]>0.0);
      indl[j]-=i;
      i=-1;
      do{
	i+=1;
	if((indr[j]+i)==(indl[j+1]+1)) break;
      }while(diffwss[indr[j]+i-1]<0.0);
      indr[j]+=i;
    }
    if(indl[j]>(indr[j-1]+1)){
      i=0;
      do{
	i+=1;
	if((indl[j]-i) <= (indr[j-1]+1)) break;
	if((xex[j-1]-x[indl[j]-i]) > halfmax) break;
      }while(((diffwss[indl[j]-i-1]-diffwss[indl[j]-i]) > 0.0) || (fabs(diffwss[indl[j]-i-1]) > *gr));
      indl[j]-=i;
    }
    if(indr[j]<(indl[j+1]-1)){
      i=1;
      do{
	i+=1;
	if((indr[j]+i) >= (indl[j+1]-1)) break;
	if(x[indr[j]+i]-xex[j-1]>halfmax) break;
      }while(((diffwss[indr[j]+i-1]-diffwss[indr[j]+i-2]) < 0.0) || (fabs(diffwss[indr[j]+i-1]) > *gr));
      indr[j]+=i;
    }
    kn=indr[j]-indl[j]+1;
    diffsort=(double *)calloc(kn,sizeof(double));
    for(k=0;k<kn;k++)
      diffsort[k]=fabs(diffwss[indl[j]-1+k]);
    qsort(diffsort, kn,sizeof(double), sortfunction);
/*     if(diffsort[0]<*gr) */
/*       exent[j]=0; */
    free(diffsort);
  }
}


void indextremw(int *nex,
		int *indl,
		int *indr,
		int *ind)
{
  int j,i,k;
  k=-1;
  for(i=0;i<*nex;i++)
    for(j=indl[i];j<=indr[i];j++){
      k+=1;
      ind[k]=j;
    }
  *nex=k;
}

void basiserg(int *np,
	      int *n,
	      double *data_x,
	      double *data_y,
	      double *basisl,
	      double *peaks,
	      int *indl,
	      int *indr,
	      int *pindl,
	      int *pindr)
{
  double x1, x2, y1, y2, a, b;
  int i,j;
  j=0;
  indl[0]=pindl[0];
  indr[0]=pindr[0];
  for(i=1;i<*np;i++){
    if(pindl[i]>indr[j]+1){
      j+=1;
      indl[j]=pindl[i];
      indr[j]=pindr[i];
    }
    else indr[j]=pindr[i];
  }
  *np=j+1;
  for(i=0;i<*np;i++){
    x1=data_x[indl[i]-2];
    x2=data_x[indr[i]];
      y1=basisl[indl[i]-2];
      y2=basisl[indr[i]];
      Rprintf("x1: %lf, x2: %lf, y1: %lf, y2: %lf\n",x1,x2,y1,y2);
      a=(y1-y2)/(x1-x2);
      b=y1-a*x1;
      for(j=(indl[i]-1);j<indr[i];j++){
	basisl[j]=a*data_x[j]+b;
	peaks[j]=data_y[j]-basisl[j];
      }
  }
}


void pkcnt(int *indl,
	   int *indr,
	   int *nind,
	   int *pkloc,
	   int *pks,
	   double *fn,
	   int *npks)
{
  int i,j,mima;
  if(fn[1]>fn[0])
    mima = 1;
  else mima = -1;
  j=0;
  for(i=0;i<*nind;i++){
    npks[i]=0;
    while(pkloc[j]<indl[i]){
      if(mima==1)
	mima=-1;
      else
	mima=1;
      j++;
    }
    while((j<pks[0])&&(pkloc[j]>=indl[i])&&(pkloc[j]<=indr[i])){
      if(mima==1){
	npks[i]++;
	mima=-1;
      }
      else mima=1;
      j++;
    }
  }
}


void multiwdwrchngd(double *y, int *n, double thresh, int *firstwidth, int lastwidth)
{
  int j, actwidth, leftind, rightind;
  double *ysum;
  
  ysum=(double *)calloc((*n+1),sizeof(double));
  ysum[0]=0;
  for(j=1;j<=*n;j++)
    ysum[j]=ysum[j-1]+y[j-1];
  for(j=0;j<*n;j++)
    y[j]=0.0;
  for(actwidth=*firstwidth;actwidth<=lastwidth;actwidth*=2) 
    for(leftind=0,rightind=actwidth;leftind<*n;leftind=rightind,rightind+=actwidth){ 
      if(rightind>*n) rightind= *n;
      if(fabs((ysum[rightind]-ysum[leftind])/sqrt((double)(rightind-leftind)))>thresh)
	for(j=leftind;j<rightind;j++)
	  y[j]=1.0;
    } 
  free(ysum);
}



void cholesky(double *A, int *n, double *L,int *FFF)
{
  double a[*n][3],l[*n][3];
  int i,j,k,p;
  a[0][0]=0.0;
  a[1][0]=0.0;
  a[0][1]=0.0;
  l[0][0]=0.0;
  l[1][0]=0.0;
  l[0][1]=0.0;
  for(k=0;k<*n;k++){
    a[k][2]=A[k*(*n)+k];
    if(k>0){
	a[k][1]=A[k*(*n)+(k-1)];
      if(k>1)
	  a[k][0]=A[k*(*n)+(k-2)];
    }
  }
  for(k=0;k<*n;k++){
    if(a[k][2]<TOL){
      Rprintf("a[%i][2]: %lf\n", k,a[k][2]);
      Rprintf("Nicht lösbar!\n");
      *FFF=1;
      break;
    }
    l[k][2]=sqrt(a[k][2]);
    if((k+2)<*n) p=k+2;
    else p=*n-1;
    for(i=k+1;i<=p;i++){
      l[i][k-i+2]=a[i][k-i+2]/l[k][2];
      for(j=k+1;j<=i;j++)
	a[i][j-i+2]=a[i][j-i+2]-l[i][k-i+2]*l[j][k-j+2];
    }
  }
  for(k=0;k<*n;k++){
    L[k*(*n)+k]=l[k][2];
    if(k<(*n-1)){
      L[(k+1)**n+k]=l[k+1][1];
      if(k<(*n-2))
	L[(k+2)**n+k]=l[k+2][0];
    }
  }
}


void vorwaerts(double *L, int *n, double *QTY, double *CCC)
{
  int i,k;
  if(fabs(L[0])<TOL)
      Rprintf("Nicht loesbar! L[0]=0.0\n");
  else{
    CCC[0]=QTY[0]/L[0];
    for(k=1;k<*n;k++)
	if(fabs(L[k**n+k])<TOL){
	    Rprintf("Nicht loesbar! L[%i][%i]=0.0\n",k,k);
	break;
      }
      else{
	CCC[k]=0.0;
	for(i=0;i<=(k-1);i++)
	  CCC[k]+=L[k**n+i]*CCC[i];
	CCC[k]=1/L[k**n+k]*(QTY[k]-CCC[k]);
      }
  }
}


void rueckwaerts(double *L, int *n, double *GGG, double *CCC)
{
  int i,k,l;
  if(fabs(L[(*n-1)**n+(*n-1)])<TOL)
    Rprintf("Nicht loesbar! L[%i]=0.0\n",(*n-1)**n+(*n-1));
  else{
    GGG[*n-1]=CCC[*n-1]/L[(*n-1)**n+(*n-1)];
    for(k=*n-2;k>=0;k=k-1)
	if(fabs(L[k**n+k])<TOL){
	    Rprintf("Nicht loesbar! L[%i][%i]=0.0\n",k,k);
	break;
      }
      else{
	GGG[k]=0.0;
	if(k-1<0) l=0;
	else l=k-1;
	for(i=l;i<=(*n-1);i++)
	  GGG[k]+=L[i**n+k]*GGG[i];
	GGG[k]=1/L[k**n+k]*(CCC[k]-GGG[k]);
      }
  }
}


void neben(double *Q, double *w, double *a, double *R,int *n,double *erg)
{
  int i,j,k,r,ll,rr;
  for (k=0;k<(*n-2);k++){
    if((k+2)>(*n-3)) r=*n-3;
    else r=k+2;
    for(j=k;j<=r;j++){
      erg[k*(*n-2)+j]=0.0;
      if(k>j) ll=k;
      else ll=j;
      if(*n>(j+2)) rr=*n;
      else 
	if(*n>(k+2)) rr=*n;
	else 
	  if((j+2)>(k+2)) rr=j+2;
	  else rr=k+2;
      for(i=ll;i<rr;i++)
	erg[k*(*n-2)+j]+=w[i]*Q[k*(*n)+i]*Q[j*(*n)+i];
      erg[k*(*n-2)+j]*=*a;
      erg[k*(*n-2)+j]+=R[j*(*n-2)+k];
      if (j>k)
	erg[j*(*n-2)+k]=erg[k*(*n-2)+j];
    }
  }
}


void wsspoisschngd(double *x, 
		   double *y, 
		   int *n, 
		   double *w, 
		   double *thresh, 
		   double *reg, 
		   int *glob, 
		   int *mr, 
		   int *nit,
		   double *q,
		   double *sqfn,
		   int *shrtint,
           int *memok)
{
  int i,j,k, firstwdth;
  double h[*n-1], *qty, *Q, *R, *T, *w1,*L, *D, *G, *res;
  double *a;
  int *Z, STP, STT, *F, *N, counter;
  a=(double *)calloc(1,sizeof(double));
  Z=(int *)calloc(1,sizeof(int));
  F=(int *)calloc(1,sizeof(int));
  N=(int *)calloc(1,sizeof(int));
  qty = (double *)calloc((*n-2),sizeof(double));
  Q = (double *)calloc(*n*(*n-2),sizeof(double));
  R = (double *)calloc((*n-2)*(*n-2),sizeof(double));
  T = (double *)calloc((*n-2)*(*n-2),sizeof(double));
  L = (double *)calloc((*n-2)*(*n-2),sizeof(double));
  w1 = (double *)calloc(*n,sizeof(double));
  D = (double *)calloc((*n-2),sizeof(double));
  G = (double *)calloc((*n-2),sizeof(double));
  res= (double *)calloc(*n,sizeof(double));

  *memok = (a==NULL)+(Z==NULL)+(F==NULL)+(N==NULL)+(qty==NULL)+(Q==NULL)+(R==NULL)+(T==NULL)+(L==NULL)+(w1==NULL)+(D==NULL)+(G==NULL)+(res==NULL);
  if (*memok!=0) {
                 Rprintf("Not enough memory for spline approximation!");
                 }

  if (*memok==0) {

  counter=0;
  *a = 0.01;
  *Z = 1;
  STP=1;
  for(i=0;i<(*n-1);i++)
    h[i] = x[i+1]-x[i];  
  for(i=0;i<(*n-2);i++)
    qty[i] = (y[i+2]-y[i+1])/h[i+1]-(y[i+1]-y[i])/h[i];
  for(i=0;i<(*n-2);i++){
    Q[i**n+i] = 1.0/h[i];
    Q[i**n+i+1] = -1.0/h[i]-1.0/h[i+1];
    Q[i**n+i+2] = 1.0/h[i+1];
    if(i<(*n-3)){
      R[i*(*n-2)+i] = 1.0/3.0*(h[i]+h[i+1]);
      R[i*(*n-2)+i+1] = R[(i+1)*(*n-2)+i] = 1.0/6.0*h[i+1];
    }
  }
  R[(*n-3)*(*n-2)+(*n-3)]=1.0/3.0*(h[*n-3]+h[*n-2]);
  *N = *n-2;
  firstwdth=1;
  if(*mr==1){
    if(*shrtint==1){
    do{
      counter+=1;
      Rprintf("Iteration: %i\n", counter); 
      *F=0;
      STP=0;
      for(i=0;i<*n;i++)
	w1[i] = 1.0/w[i];
      neben(Q,w1,a,R,n,T);
      cholesky(T, N, L, F);
      if(*F==1) break;
      vorwaerts(L,N, qty,D);
      rueckwaerts(L, N, G, D);
      for(j=2;j<(*n-2);j++){
	reg[j] = 0.0;
	for(i=j-2;i<=j;i++){
	  reg[j]+=Q[i**n+j]*G[i];
	}
	reg[j] = y[j]-w1[j]**a*reg[j];
      }
      reg[0] = y[0]-w1[0]**a*Q[0]*G[0];
      reg[1] = y[1]-w1[1]**a*(Q[1]*G[0]+Q[1**n+1]*G[1]);
      reg[*n-2] = y[*n-2]-w1[*n-2]**a*(Q[(*n-4)**n+*n-2]*G[*n-4]+Q[(*n-3)**n+*n-2]*G[*n-3]);
      reg[*n-1] = y[*n-1]-w1[*n-1]**a*Q[(*n-3)**n+*n-1]*G[*n-3];

       for(i=0;i<*n;i++) 	res[i] = (y[i]-reg[i])/sqfn[i];
       multiwdwrchngd(res,n,*thresh, Z,32);
       STT=0;
       if(*glob>0)
	for(i=0;i<*n;i++){
	  if(res[i]==1.0){
	    for(k=0;k<*n;k++)
	      w[k]*=*q;
	    STP=1;
	    STT=1;
	  }
	  if(STT==1)
	    break;
	}
      else
	for(i=0;i<*n;i++){
	  if(res[i]==1.0){
	    w[i]*=*q;
	    STP=1;
	  }
	}
    }while (STP==1);
    }
    do{
      counter+=1;
      Rprintf("Iteration: %i\n", counter); 
      *F=0;
      STP=0;
      for(i=0;i<*n;i++)
	w1[i] = 1.0/w[i];
      neben(Q,w1,a,R,n,T);
      cholesky(T, N, L, F);
      if(*F==1) break;
      vorwaerts(L,N, qty,D);
      rueckwaerts(L, N, G, D);
      for(j=2;j<(*n-2);j++){
	reg[j] = 0.0;
	for(i=j-2;i<=j;i++){
	  reg[j]+=Q[i**n+j]*G[i];
	}
	reg[j] = y[j]-w1[j]**a*reg[j];
      }
      reg[0] = y[0]-w1[0]**a*Q[0]*G[0];
      reg[1] = y[1]-w1[1]**a*(Q[1]*G[0]+Q[1**n+1]*G[1]);
      reg[*n-2] = y[*n-2]-w1[*n-2]**a*(Q[(*n-4)**n+*n-2]*G[*n-4]+Q[(*n-3)**n+*n-2]*G[*n-3]);
      reg[*n-1] = y[*n-1]-w1[*n-1]**a*Q[(*n-3)**n+*n-1]*G[*n-3];
      
      for(i=0;i<*n;i++) 	res[i] = (y[i]-reg[i])/sqfn[i];
      
      multiwdwrchngd(res,n,*thresh, Z,*n);
      STT=0;
      if(*glob>0)
	for(i=0;i<*n;i++){
	  if(res[i]==1.0){
	    for(k=0;k<*n;k++)
	      w[k]*=*q;
	    STP=1;
	    STT=1;
	  }
	  if(STT==1)
	    break;
	}
      else
	for(i=0;i<*n;i++){
	  if(res[i]==1.0){
	    w[i]*=*q;
	    STP=1;
	  }
	}
    }while (STP==1);
  }
  else
    do{
      counter+=1;
      Rprintf("Iteration: %i\n", counter); 
      *F=0;
      STP=0;
      for(i=0;i<*n;i++)
	w1[i] = 1.0/w[i];
      neben(Q,w1,a,R,n,T);
      cholesky(T, N, L, F);
      if(*F==1) break;
      vorwaerts(L,N, qty,D);
      rueckwaerts(L, N, G, D);
      for(j=2;j<(*n-2);j++){
	reg[j] = 0.0;
	for(i=j-2;i<=j;i++){
	  reg[j]+=Q[i**n+j]*G[i];
	}
	reg[j] = y[j]-w1[j]**a*reg[j];
      }
      reg[0] = y[0]-w1[0]**a*Q[0]*G[0];
      reg[1] = y[1]-w1[1]**a*(Q[1]*G[0]+Q[1**n+1]*G[1]);
      reg[*n-2] = y[*n-2]-w1[*n-2]**a*(Q[(*n-4)**n+*n-2]*G[*n-4]+Q[(*n-3)**n+*n-2]*G[*n-3]);
      reg[*n-1] = y[*n-1]-w1[*n-1]**a*Q[(*n-3)**n+*n-1]*G[*n-3];

      /* HIER POISSON-AENDERUNG!!!!! */
      for(i=0;i<*n;i++) 	res[i] = (y[i]-reg[i])/sqfn[i];

      multiwdwrchngd(res,n,*thresh, Z,*n);
      STT=0;
      if(*glob>0){
	if(counter<(*nit-1))
	  for(i=0;i<*n;i++){
	    if(res[i]==1.0){
	      for(k=0;k<*n;k++)
		w[k]*=*q;
	      STP=1;
	      STT=1;
	    }
	    if(STT==1)
	      break;
	  }
	}
	else
	if(counter<(*nit-1))
	for(i=0;i<*n;i++){
	  if(res[i]==1.0){
	    w[i]*=*q;
	    STP=1;
	  }
	}
    }while(counter<*nit);

    }                      

  free(N);
  free(a);
  free(Z);
  free(F);
  free(qty);
  free(Q);
  free(R);
  free(T);
  free(w1);
  free(L);
  free(D);
  free(G);
  free(res);
}
