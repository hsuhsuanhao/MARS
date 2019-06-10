#ifndef OPT_h
#define OPT_h
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "MATRIX.h"
#include "ATOM.h"
#include <fstream>
#include <string>
using namespace std;
class OPT {
	public :
		OPT (){
			data.set_up();
			Etot=0.0;
			act=0;
			tot=0;
			count=0;
			id=NULL;
			table=NULL;
		}
		int optimize();
		int bnd(int, int);  // energy of bond
		int ang(int , int ,int); // energy of angle
		int validate();
		void output(ostream &);
		double Etot;
		int act,tot;
		Matrix x0,f,df;
		POOL data;
		int **table, *id, num,count;
};
#endif
