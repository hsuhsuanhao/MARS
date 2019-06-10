#include "ATOM.h"
#include "MOLECULE.h"
#include "MATRIX.h"
#include "OPT.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>

using namespace std;
void mark(ostream &);

int main(int argc, char **argv) {
	mark(cout);

	int i,j,k,dd3=0,enc=0;
	if (argc==1) {
		cout<<"Error : no input file. Please see README.txt"<<endl;
		exit(1);
	}

	string infile=argv[1];
	vector<string> command;
	if (argc==3) command.push_back(argv[2]);
	else if (argc==4) {
		command.push_back(argv[2]);
		command.push_back(argv[3]);
	}
	for (i=0;i<command.size();i++) {
		if (command.at(i)=="gen3d"||command.at(i)=="gen3D"||command.at(i)=="Gen3d"||command.at(i)=="Gen3D") dd3=1;
		else if (command.at(i)=="genmds"||command.at(i)=="GenMDS" ||command.at(i)=="genMDS") enc=1;
	}

	ifstream inf;
	MOLECULE *gs;
	gs = new MOLECULE [2];
	string smiles;
	inf.open(infile.c_str());
	ofstream log("molecule.log");
	mark(log);

	i=0;
	k=1;
	cout<<"Reading molecules from input file "<<infile.c_str()<<endl;
	while (!inf.eof()) {
		smiles="";
		inf>>smiles>>ws;
		if (smiles=="" && i==1) {
			cout<<"Warning : for combination and crossover, you should obtain 2 molecules."<<endl;
			log<<"Warning : for combination and crossover, you should obtain 2 molecules."<<endl;
			k=0;
			break;
		}
		else if (smiles=="" && i==0) {
			log<<"Error :  "<<infile<<" is empty."<<endl;
			cout<<"Error : "<<infile<<" is empty."<<endl;
			exit(1);
		}
		gs[i].smiles=smiles;
		gs[i].input();
		i++;
	}
	cout<<"Done. "<<i<<" molecules read."<<endl;
	inf.close();

	ofstream out;
	if (k) {
		cout<<"Performing combination operation on "<<gs[0].smiles<<" and "<<gs[1].smiles<<endl;
		out.open("combination.txt");
		mark(out);
		for (i=0;i<gs[0].Cindex.size();i++) {
			for (j=0;j<(gs[1].Cindex.size());j++) {
				if (gs[0].combine(gs[1],i,j)) {
					out<<gs[0].molesmi<<endl;
					if (enc) {
						gs[0].outmds(out);
					}
					if (dd3) {
						gs[0].mds23d(out);
					}
					gs[0].input();
				}
			}
		}
		out.close();
		cout<<"Performing crossover operation on "<<gs[0].smiles<<" and "<<gs[1].smiles<<endl;
		out.open("crossover.txt");
		mark(out);
		for (i=1;i<gs[0].Cindex.size();i++) {
			for (j=1;j<gs[1].Cindex.size();j++) {
				if (gs[0].crossover(gs[1],i,j)) {
					gs[0].mds2smi();
					gs[1].mds2smi();
					out<<gs[0].molesmi<<endl;
					if (enc) {
						gs[0].outmds(out);
					}
					if (dd3) gs[0].mds23d(out);
					out<<gs[1].molesmi<<endl;
					if (enc) {
						gs[1].outmds(out);
					}
					if (dd3) gs[1].mds23d(out);
					gs[0].input();
					gs[1].input();
				}
			}
		}
		out.close();
	}

	cout<<"Performing subtraction operation on "<<gs[0].smiles<<endl;
	out.open("subtraction.txt");
	mark(out);
	for (i=1;i<gs[0].Cindex.size();i++) {
		gs[0].subtract(i);
		gs[0].mds2smi();
		out<<gs[0].molesmi<<endl;
		if (enc) {
			gs[0].outmds(out);
		}
		if (dd3) gs[0].mds23d(out);
		gs[0].input();
	}
	out.close();

	cout<<"Performing addition operation on "<<gs[0].smiles<<endl;
	out.open("addition.txt");
	mark(out);
	for (i=0;i<gs[0].Cindex.size();i++) {
		for (j=1;j<36;j++) {
			if (gs[0].add(i,j)) {
				gs[0].mds2smi();
				out<<gs[0].molesmi<<endl;
				if (enc) {
					gs[0].outmds(out);
				}
				if (dd3) gs[0].mds23d(out);
				gs[0].input();
			}
		}
	}
	out.close();

	cout<<"Performing exchange operation on "<<gs[0].smiles<<endl;
	out.open("exchange.txt");
	mark(out);
	int m=0;
	for (i=1;i<gs[0].Cindex.size();i++) {
		for (j=1;j<36;j++) {
			for (k=-2;k<3;k++) {
				for (m=1;m<36;m++) {
					if (gs[0].exchange(i,j,m,k)) {
						gs[0].mds2smi();
						out<<gs[0].molesmi<<endl;
						if (enc) gs[0].outmds(out);
						if (dd3) gs[0].mds23d(out);
						gs[0].input();
					}
				}
			}
		}
	}
	out.close();

	cout<<"Performing ring operation on "<<gs[0].smiles<<endl;
	out.open("ring.txt");
	mark(out);
	for (i=0;i<gs[0].Cindex.size();i++) {
		for (j=i+1;j<gs[0].Cindex.size();j++) {
			if (gs[0].ring(i,j)) {
				gs[0].mds2smi();
				out<<gs[0].molesmi<<endl;
				gs[0].outmds(out);
				if (dd3) gs[0].mds23d(out);
				gs[0].input();
			}
		}
	}

	cout<<"Creating benzene, c1ccccc1 from methane, C"<<endl;
	out.close();
	gs[0].smiles="methane";
	gs[0].input();
	out.open("benzene.txt");
	mark(out);
	out<<gs[0].molesmi<<endl;
	if (dd3) gs[0].mds23d(out);
	gs[0].add(0,1);
	gs[0].mds2smi();
	out<<gs[0].molesmi<<endl;
	if (dd3) gs[0].mds23d(out);
	gs[0].exchange(1,2,2,1);
	gs[0].mds2smi();
	out<<gs[0].molesmi<<endl;
	if (dd3) gs[0].mds23d(out);
	gs[0].add(1,2);
	gs[0].mds2smi();
	out<<gs[0].molesmi<<endl;
	if (dd3) gs[0].mds23d(out);
	gs[0].add(2,2);
	gs[0].mds2smi();
	out<<gs[0].molesmi<<endl;
	if (dd3) gs[0].mds23d(out);
	gs[0].add(3,2);
	gs[0].mds2smi();
	out<<gs[0].molesmi<<endl;
	if (dd3) gs[0].mds23d(out);
	gs[0].add(4,2);
	gs[0].mds2smi();
	out<<gs[0].molesmi<<endl;
	if (dd3) gs[0].mds23d(out);
	gs[0].ring(0,5);
	gs[0].mds2smi();
	out<<gs[0].molesmi<<endl;
	if (dd3) gs[0].mds23d(out);
	out.close();
	cout<<"Job completed!"<<endl;
	cout<<"Results can be found in corresponding .txt files. (combination.txt, crossover.txt, addition.txt, subtraction.txt, exchange, ring.txt, and benzene.txt)"<<endl;
	delete [] gs;
	return 1;
}

void mark(ostream& out) {
	    out <<"***************************************************************************************"<<endl
		<<"                   MARS: Molecular Assembling and Representation Suite                 "<<endl
		<<"                                      Hsuan Hao Hsu                                    "<<endl
		<<"                                    Chen Hsuan Huang                                   "<<endl
		<<"                            Shiang Tai Lin (stlin@ntu.edu.tw)                          "<<endl
		<<"                      Computational Molecular Engineering Laboratory                   "<<endl
		<<"      Department of Chemical Engineering, National Taiwan University, Taipei, Taiwan   "<<endl
		<<"                                   Copyright (c) 2019                                  "<<endl
		<<"                                   All rights Reserved                                 "<<endl
		<<"***************************************************************************************"<<endl;
	return;
}
