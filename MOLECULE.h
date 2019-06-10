/***************************************************************************************
                          MARS for Computer-aided Molecular Design
                                      Hsuan Hao Hsu
                                    Chen Hsuan Huang
                            Shiang Tai Lin (stlin@ntu.edu.tw)
      Department of Chemical Engineering, National Taiwan University, Taipei, Taiwan
                                   Copyright (c) 2019
                                   All rights Reserved
***************************************************************************************/
#ifndef MOLECULE_h
#define MOLECULE_h
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
#include "ATOM.h"
#include "MATRIX.h"
#include "OPT.h"
using namespace std;
class MOLECULE
{
	public:
		MOLECULE() {
			data.set_up();
			molesmi="";
			parent=-1;
			if_circle=0;
			dist=NULL;
			connect=NULL;
			order=NULL;
			atm=NULL;
			natom=0;
			Pindex.resize(0);
			Cindex.resize(0);
			Mindex.resize(0);
			Cyindex.resize(0);
			Rindex.resize(0);
			for (int i=0;i<1024;i++) {
				Bindex[i].resize(0);
				atomsmi[i].resize(0);
			}
		}

		~MOLECULE() {
				int i;
				vector<int>().swap(Pindex);
				vector<int>().swap(Cindex);
				vector<int>().swap(Mindex);
				vector<int>().swap(Cyindex);
				vector<int>().swap(Rindex);
				for (i=0;i<3;i++) vector<double>().swap(xyz[i]);
				for (i=0;i<1024;i++) {
					vector<int>().swap(Bindex[i]);
					vector<int>().swap(atomsmi[i]);
				}
		}
		POOL data;
		int chg;
		void mds23d(ostream &outs=cout);
		int ring(int,int);
		void mds2smi();
		void subtract(int);
		void neutralize();
		int add(int,int);
		void reset();
		void clear();
		void empty();
		void clean();
		void wipe();
		void outmds(ostream &outs=cout);
		int exchange(int,int,int,int);
		void read(string);
		void smi2gjf();
		void init();
		void init(int);
		void cal_r();
		void input();
		void check_bnd();
		void check_bnd(int);
		void smi2cod();
		void report();
		int rd_nps(int);
		int crossover(MOLECULE &, int, int);
		int combine(MOLECULE &, int, int);						// give a exact point to combine
		void recode(int);
		int chk_chg(MOLECULE &);
		int rechg();
		double prob();
		vector<int> Pindex; 															// parent index
		vector<int> Cindex; 															// child index
		vector<int> Mindex; 															// molecule index
		vector<int> Bindex[1024]; 														// to measure bond to 0, do not output.
		vector<int> Rindex;
		vector<int> Cyindex;
		vector<int> atomsmi[1024];
		string smiles;
		int natom;
		double **dist;
		int **connect;
		int **order;
		DEATOM *atm;
		int parent;
		int if_circle;
		double p_sel;
		vector <double> xyz[3]; 	
		string molesmi; 								// smiles representation
		string ionic;
	private:
};
#endif
