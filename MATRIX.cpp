#include "MATRIX.h"
#include <cmath>
using namespace std;

Matrix::Matrix()
: nRows(0), nCols(0),Mtrx(0), rank(0) ,det(0)
{ }

Matrix::Matrix( int rows, int cols )
: nRows(rows), nCols(cols), Mtrx(rows)
{
   int r,c;
   for (r=0; r<nRows; r++)
      Mtrx[r].resize(cols);
}

Matrix::~Matrix() {
	vector<vector<double> >().swap(Mtrx);
}
ostream& operator<< (ostream& outs, Matrix m2)
{
    m2.display(outs);
    return outs;
}

int Matrix::row()
{
   return nRows;
}

int Matrix::col()
{
   return nCols;
}

double Matrix::dtm()
{
   inverse();
   return det;
}


int Matrix::assign(int r,int c,double val)
{
   Mtrx[r][c]=val;
   return 0;
}

int Matrix::resize(int newRows, int newCols)
{
	Mtrx.resize(newRows);
    for(int r=0; r < newRows; r++)
        Mtrx[r].resize(newCols);
    nRows = newRows;
    nCols = newCols;
	int i,j;
	for (i=0;i<nRows;i++) {
		for (j=0;j<nCols;j++) {
			Mtrx[i][j]=0.00000;
		}
	}
    return 0;
}

int Matrix::display(ostream& outs)
{
    int i,j;
    outs<<"This is a "<<nRows<<"x"<<nCols<<" matrix."<<endl;
    for(i=0;i<nRows;i++) {
        for(j=0;j<nCols;j++) {
            outs.width(9); outs.precision(4); outs.flags(ios::fixed);
            outs<<Mtrx[i][j]<<" ";
        }
        outs<<endl;
    }
    outs<<endl;
    return 0;
}

vector<double>&Matrix::operator[] (int row)
{
   check_range(row);
   return Mtrx[row];
}

int Matrix::check_range (int row)
{
   if (row< 0 || row>= nRows ) {
      cerr << "\n***ERROR: index [" << row
           << "] out of range (" << 0
           << ".." << (nRows-1) << ")\n";
      }
   return 0;
}


int Matrix::ludcmp(Matrix &a,int *indx,double &d)
{
    int n=rank;
    int i,imax,j,k;
    double big,dum,sum,temp;
    double *vv=new double[n];
    double TINY=1.0e-20;

    d=1.0;
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++) {
            if ( ( temp = fabs(a[i][j]) ) > big ) big=temp;
        }
        if (big == 0.0) cerr<<"Singular matrix in routine ludcmp";
        vv[i]=1.0/big;
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            sum=a[i][j];
            for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<n;i++) {
                sum=a[i][j];
                for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
                a[i][j]=sum;
                if ( (dum=vv[i]*fabs(sum)) >= big) {
                    big=dum;
                    imax=i;
                }
        }
        if (j != imax) {
                for (k=0;k<n;k++) {
                    dum=a[imax][k];
                    a[imax][k]=a[j][k];
                    a[j][k]=dum;
                }
                d = -d;
                vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<n;i++) a[i][j] *= dum;
        }
    }
    delete [] vv;
    return 0;
}

int Matrix::lubksb(Matrix &a,int *indx,double *b)
{
    int n=rank;
    int i,ii=0,ip,j;
    double sum;

    for (i=0;i<n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii+1) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
    }
    for (i=n-1;i>=0;i--) {
        sum=b[i];
        for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
    return 0;
}
Matrix Matrix::inverse()
{
   Matrix mInv;
   if(nRows!=nCols) {
      cout<<"error: not a square matrix, inversion not possible"<<endl;
   } else {
      rank=nRows;
      int *indx,i,j;
      double *col;
      Matrix Ac;
      indx=new int[rank];
      col=new double[rank];

      Ac=*this;

      ludcmp(Ac,indx,det);

      for(i=0;i<rank;i++) det *= Ac[i][i];

      mInv.resize(nRows,nCols);
      for(j=0;j<rank;j++) {
        for(i=0;i<rank;i++) col[i]=0.0;
           col[j]=1.0;
           lubksb(Ac,indx,col);
           for(i=0;i<rank;i++) mInv[i][j]=col[i];
         }
       }

   return mInv;
}
Matrix Matrix::operator+ (Matrix m2)
{
   Matrix sum(nRows,nCols);
   if(nRows!=m2.nRows || nCols!=m2.nCols) {
       cout<<"addition not possible"<<endl;
   } else {
       int i,j;
       for(i=0;i<nRows;i++) {
          for(j=0;j<nCols;j++)  sum[i][j]=Mtrx[i][j]+m2[i][j];
       }
   }
   return sum;
}

Matrix Matrix::operator- (Matrix m2)
{
   Matrix sum(nRows,nCols);
   if(nRows!=m2.nRows || nCols!=m2.nCols) {
       cout<<"addition not possible"<<endl;
   } else {
       int i,j;
       for(i=0;i<nRows;i++) {
          for(j=0;j<nCols;j++)  sum[i][j]=Mtrx[i][j]-m2[i][j];
       }
   }
   return sum;
}

Matrix Matrix::operator* (Matrix m2)
{
   Matrix sum(nRows,m2.nCols);
   if(nCols!=m2.nRows) {
       cout<<"multiplication not possible"<<endl;
   } else {
       int i,j,k;
       for(i=0;i<nRows;i++) {
           for(j=0;j<m2.nCols;j++) {
               sum[i][j]=0;
               for(k=0;k<nCols;k++) sum[i][j]+= Mtrx[i][k]*m2[k][j];
           }
       }
   }
   return sum;
}

Matrix Matrix::operator* (double a) {
	Matrix sum(nRows,nCols);
	int i,j;
	for (i=0;i<nRows;i++) {
		for (j=0;j<nCols;j++) {
			sum[i][j]=Mtrx[i][j]*a;
		}
	}
	return sum;
}
Matrix Matrix::operator^ (int p)
{
   Matrix sum;
   if(nCols!=nRows) {
       cout<<"power not possible"<<endl;
   } else if (p>0){
       sum=*this;
       int i;
       for(i=1;i<p;i++) sum=sum*(*this);
   } else if(p==-1) {
       sum=inverse();
   } else {
      cout<<"power of "<<p<<"not allowed"<<endl;
   }
   return sum;
}

void Matrix::rotate(Matrix &a, Matrix &b, double deg) {
	int i,j,k;
	double x,y;
	
	x=sqrt(pow(a[0][0],2)+pow(a[1][0],2)+pow(a[2][0],2));
	y=sqrt(pow(a[1][0],2)+pow(a[2][0],2));
	Matrix Rx(3,3),Ry(3,3),Rz(3,3),rx(3,3),ry(3,3);
	for (i=0;i<3;i++) {
		for (j=0;j<3;j++) {
			Rx[i][j]=0.0;
			Ry[i][j]=0.0;
			Rz[i][j]=0.0;
			rx[i][j]=0.0;
			ry[i][j]=0.0;
		}
	}
	Rx[0][0]=rx[0][0]=1.0;
	Rx[1][1]=Rx[2][2]=rx[1][1]=rx[2][2]=a[2][0]/y;
	Rx[1][2]=rx[2][1]=-a[1][0]/y;
	Rx[2][1]=rx[1][2]=a[1][0]/y;
	Ry[0][0]=Ry[2][2]=ry[0][0]=ry[2][2]=y/x;
	Ry[1][1]=ry[1][1]=1.0;
	Ry[0][2]=ry[2][0]=-a[0][0]/x;
	Ry[2][0]=ry[0][2]=a[0][0]/x;
	Rz[0][0]=Rz[1][1]=cos(deg);
	Rz[0][1]=-sin(deg);
	Rz[1][0]=sin(deg);
	Rz[2][2]=1.0;
	b=rx*ry*Rz*Ry*Rx*b;
	return;
}

int Matrix::sol(Matrix &ans,Matrix &nor,Matrix &x0) {
	double err;
	int i,j;
	Matrix x1(ans.nRows,ans.nCols),f(ans.nRows,ans.nCols),df(ans.nRows,ans.nRows);
	while (1) {
		func(f,nor,x0,ans);
		dfun(df,nor,x0,ans);
		x1=ans-(df^-1)*f;
		err=0.0;
		for (i=0;i<x1.nRows;i++) {
			for (j=0;j<x1.nCols;j++) {
				err+=fabs(x1[i][j]-ans[i][j]);
			}
		}
		ans=x1;
		if (err<1E-5) break;
	}
	return 1;
}

int Matrix::func(Matrix &f, Matrix &n, Matrix &x0, Matrix &x1) {
    int i;
    double sum=0;
    for (i=0;i<3;i++) sum+=(x1[i][0]-x0[i][0])*(x1[i][0]-x0[i][0]);
    sum-=1.0;
	f[0][0]=n[0][0]*(x1[0][0]-x0[0][0])+n[1][0]*(x1[1][0]-x0[1][0])+n[2][0]*(x1[2][0]-x0[2][0]);
	f[1][0]=sum;
	sum=0.0;
	for (i=0;i<3;i++) sum+=(x1[i][0]-x0[i][0])*n[i][0];
	f[2][0]=sum;
	return 1;
}

int Matrix::dfun(Matrix &df, Matrix &n, Matrix &x0, Matrix &x1) {
	int i,j,k;
	double dh=1E-5;
	Matrix tmp, tmpp;
	for (j=0;j<x1.nCols;j++) {
		tmp=x1;
		tmpp=x1;
		for (i=0;i<tmp.nRows;i++) {tmp[i][j]+=dh;tmpp[i][j]-=dh;}
		func(tmp,n,x0,tmp);
		func(tmpp,n,x0,tmpp);
		for (i=0;i<tmp.nRows;i++) {df[i][j]=(tmp[i][j]-tmpp[i][j])/2.0/dh;}
	}
	return 1;
}
	
