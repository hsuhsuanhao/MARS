#ifndef ATOM_h
#define ATOM_h
#include <sys/types.h>
#include <dirent.h>
#include <cerrno>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <vector>
#include <ctime>
#include <cstdio> 
#include <cctype>
#include <sstream>
#include "MATRIX.h"
using namespace std;
class ATOM
{
	public:
		ATOM() {
			order.resize(0);
		}
		~ATOM() {
			vector<int>().swap(order);
		}
		int id;
		string name,atm;
		vector<int> order;
		int nh;
		int bd[3];
		int type;
		int index;
		double rb;
		double ang0;
		int norder;
		int chg;
		int nbond; // how many char for name
};

class DEATOM {
	public :
		DEATOM() {
		}
		~DEATOM() {}
		string name;
		double r_bnd;
		int nbnd;
		double x[3];
		void find_r();
}; 

class POOL
{
	public:
		POOL(){ 
			//a=NULL;
			a.resize(0);
			if (0) set_up();
		}
		~POOL(){ 
			if (0) {
				/*
				if(a!=NULL) {
					delete [] a;
					a=NULL;
				}
				*/
				//a.clear();
				vector<ATOM>().swap(a);
			}
		}
		int num;
		//ATOM *a; //  types of atoms
		vector<ATOM> a;
		void set_up();

	private:
};
#endif
