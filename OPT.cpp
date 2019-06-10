#include "OPT.h"
#include "MATRIX.h"
#include "ATOM.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;
int OPT::optimize() {
	int i,j,k,m,n=0;
	vector<int> tmp;
	double sum=0.0;
	double tt;
	Matrix x1(num*3,1),delta(num*3,1),step(num*3,1);
	if (count) count=10000;
	else count=100;
	while (1) {
		sum=0.0;
		for (i=0;i<3*num;i++) step[i][0]=1E-2;
		for (i=0;i<3*num;i++) df[0][i]=0.0;
		Etot=0.0;
		n++;
		if (1) {
		sum=Etot=0.0;
		for (i=0;i<num;i++) {
			sum+=f[i][0];
		}
		k=0;
		tt=0.0;
    	for (i=0;i<num;i++) f[i][0]=0.0;
   			for (i=0;i<num;i++) {
        		for (j=i+1;j<num;j++) {
            		if (table[i][j]) {
                	bnd(i,j);
            	}
        	}
    	}
    	for (i=0;i<num;i++) {
        	for (j=0;j<num;j++) {
            	if (table[i][j] && i!=j) tmp.push_back(j);
        	}
        	if (tmp.size()>1) {
            	for (k=0;k<tmp.size();k++) {
                	for (m=k+1;m<tmp.size();m++) {
                        ang(tmp.at(k),i,tmp.at(m));
                    }
                }
        	}
        	tmp.clear();
   		}
		for (i=0;i<num;i++) Etot+=f[i][0];
		for (i=0;i<3*num;i++) {
			if (fabs(df[0][i])>tt) tt=fabs(df[0][i]);
		}
		if (1) {
			for (j=0;j<3*num;j++) {
				step[j][0]/=(tt/0.5);
			}
		}
		for (i=0;i<3*num;i++) {
			delta[i][0]=step[i][0]*df[0][i];
		}
		x1=x0-delta;
		n++;
		x0=x1;
		for (i=0;i<num;i++) sum+=f[i][0];
		}
		if (fabs(sum-Etot) <1E-5 && n>1000) break;
		else if (n>count) break;
	}
	return 1;
}



int OPT::ang(int l,int m,int n) {
	int i,j,k;
	double ang0=data.a[id[m]].ang0;
	double dot=0.0,deg;
	double r1[3],r2[3],R1=0.0,R2=0.0,C=120.0,f1[3],f2[3];
	for (i=0;i<3;i++) {
		r1[i]=x0[l+num*i][0]-x0[m+num*i][0];
		r2[i]=x0[n+num*i][0]-x0[m+num*i][0];
	}
	
	for (i=0;i<3;i++) {
		dot+=r1[i]*r2[i];
		R1+=r1[i]*r1[i];
		R2+=r2[i]*r2[i];
	}
	R1=sqrt(R1);
	R2=sqrt(R2);
	deg=dot/R1/R2; // cos theta
	for (i=0;i<3;i++) {
		f[l][0]+=C/sin(ang0)/sin(ang0)*(deg-cos(ang0))*(deg-cos(ang0));
		f[m][0]+=C/sin(ang0)/sin(ang0)*(deg-cos(ang0))*(deg-cos(ang0));
		f[n][0]+=C/sin(ang0)/sin(ang0)*(deg-cos(ang0))*(deg-cos(ang0));
	}
	
	for (i=0;i<3;i++) {
		f1[i]=C/sin(ang0)/sin(ang0)*2.0*(deg-cos(ang0))*(r2[i]/R1/R2-r1[i]*dot/R2/R1/R1/R1);
		f2[i]=C/sin(ang0)/sin(ang0)*2.0*(deg-cos(ang0))*(r1[i]/R1/R2-r2[i]*dot/R1/R2/R2/R2);
	}
	//cout<<deg<<" "<<ang0<<" "<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<endl;
	for (i=0;i<3;i++) {
		df[0][l+i*num]+=f1[i];
		df[0][n+i*num]+=f2[i];
		df[0][m+i*num]-=(f1[i]+f2[i]);
	}
	return 1;
}

int OPT::bnd(int n, int m) {
	int i,j,k;
	double r[3],R=0.0,r0=(data.a[id[n]].rb+data.a[id[m]].rb)*1.00,C=500.0;
	for (i=0;i<3;i++) r[i]=x0[m+i*num][0]-x0[n+i*num][0];
	for (i=0;i<3;i++) R+=r[i]*r[i];
	R=sqrt(R);
	for (i=0;i<3;i++) {
		f[n][0]+=0.5*C*(R-r0)*(R-r0);
		f[m][0]+=0.5*C*(R-r0)*(R-r0);
	}
	for (i=0;i<3;i++) {
		df[0][n+i*num]+=-C*(R-r0)*r[i]/R;
        df[0][m+i*num]-=-C*(R-r0)*r[i]/R;
	}
	return 1;
}

int OPT::validate() {
	int i,j,k;
	cout<<"TRUE ANSWER : "<<endl;
	for (i=0;i<num;i++) {
		for (j=0;j<num;j++) {
			cout<<table[i][j]<<" ";
		}
		cout<<endl;
	}

	cout<<"RESULT OF OPTIMIZATION :"<<endl;
	int tmp[num][num];
	double r0,r;
	for (i=0;i<num;i++) {
		for (j=0;j<num;j++) tmp[i][j]=0;
	}
	for (i=0;i<num;i++) {
		for (j=i+1;j<num;j++) {
			r0=(data.a[id[i]].rb+data.a[id[j]].rb)*1.15;
			r=0.0;
			for (k=0;k<3;k++) r+=(x0[i+k*num][0]-x0[j+k*num][0])*(x0[i+k*num][0]-x0[j+k*num][0]);
			r=sqrt(r);
			if (r<r0) tmp[i][j]=tmp[j][i]=1;
		}
	}
	for (i=0;i<num;i++) {
		for (j=0;j<num;j++) cout<<tmp[i][j]<<" ";
		cout<<endl;
	}
	k=0;
	for (i=0;i<num;i++) {
		for (j=0;j<num;j++) {
			if (tmp[i][j]==table[i][j]) {
				k=0;
			}
			else {
				k=1;
				break;
			}
		}
	}
	ofstream out("validate.log",ios::app);
	if (k) out<<1<<endl;
	else out<<0<<endl;
	out.close();
	return 1;
}		
	
void OPT::output(ostream &out) {
	int i;
	for (i=0;i<num;i++) out<<data.a[id[i]].atm<<" "<<x0[i][0]<<" "<<x0[i+num][0]<<" "<<x0[i+2*num][0]<<endl;
	return;
}
