#include "MOLECULE.h"
#include <cstdlib>
#define p_circle 0.2
//PARAMETER par;

double MOLECULE::prob() {
	double p;
	p=(double)rand()/RAND_MAX;
	return p;
}

int MOLECULE::crossover(MOLECULE &aaa,int pp, int jj) {
    int n;
    int m;
    int k=0;
    int a;
    int b;
    int y;
    int x;
    int c;
    int ref;
    int delta;
    int pos_i;
    int pos_j;
    vector<int> t,f;
    vector<int> ci_ref, mi_ref,pi_ref,ri_ref,cyi,pri;
    vector<int> cj_ref, mj_ref,pj_ref,rj_ref,cyj,prj;
    k=0;
    if (aaa.Rindex.at(jj)==Rindex.at(pp)) {
        k=1;
        m=jj;
        n=pp;
        b=aaa.Pindex[m];
        a=Pindex[n];
                y=aaa.Cindex[m];
                x=Cindex[n];
                ri_ref.push_back(Rindex[n]);
                rj_ref.push_back(aaa.Rindex[m]);
                ci_ref.push_back(x);
                cj_ref.push_back(y);
                pi_ref.push_back(a);
                pj_ref.push_back(b);
                mi_ref.push_back(Mindex[n]);
                cyi.push_back(Cyindex.at(n));
                cyj.push_back(aaa.Cyindex.at(m));
                mj_ref.push_back(aaa.Mindex[m]);
    }

    if (k==0) return 0;
    if (k==1) {
        for (n=1;n<Cindex.size();n++) {
            if (Cindex[n] > ci_ref[0]) {
                ref = ci_ref.size();
                for (m=0;m<ref;m++) {
                    if (Pindex[n] == ci_ref[m]) {
                        ci_ref.push_back(Cindex.at(n));
                        mi_ref.push_back(Mindex.at(n));
                        pi_ref.push_back(Pindex.at(n));
                        ri_ref.push_back(Rindex.at(n));
                        cyi.push_back(Cyindex.at(n));
                    }
                }
            }
        }
        for (m=1;m<aaa.Cindex.size();m++) {
            if (aaa.Cindex[m] > cj_ref[0] ) {
                ref = cj_ref.size();
                for (n=0;n<ref;n++) {
                    if (aaa.Pindex[m] == cj_ref[n]) {
                        cj_ref.push_back(aaa.Cindex.at(m));
                        mj_ref.push_back(aaa.Mindex.at(m));
                        pj_ref.push_back(aaa.Pindex.at(m));
                        rj_ref.push_back(aaa.Rindex.at(m));
                        cyj.push_back(aaa.Cyindex.at(m));
                    }
                }
            }
        }
		for (n=0;n<pri.size();n++) if (pri.at(n)) return 0;
		for (n=0;n<prj.size();n++) if (prj.at(n)) return 0;
		a=0;
		for (n=0;n<cyi.size();n++) {
			if (cyi.at(n)) {cyi[n]+=if_circle;a++;}
		}
		if (a%2 && if_circle) {
			for (n=0;n<cyi.size();n++) if (cyi.at(n)) cyi.at(n)=if_circle;
		}
		else if (a%2) {
			for (n=0;n<cyi.size();n++) cyi.at(n)=0;
		}
		a=0;
		for (n=0;n<cyj.size();n++) {
			if (cyj.at(n)) {cyj[n]+=aaa.if_circle;a++;}
		}
        if (a%2 && aaa.if_circle) {
            for (n=0;n<cyj.size();n++) if (cyj.at(n)) cyj.at(n)=aaa.if_circle;
        }
		if (a%2) {
			for (n=0;n<cyj.size();n++) cyj.at(n)=0;
		}
		
        x=ci_ref[0];
        y=cj_ref[0];
        a=pi_ref[0];
        b=pj_ref[0];
        for (n=0;n<ci_ref.size();n++) {
            for (m=0;m<Cindex.size();m++) {
                if (Cindex[m] == ci_ref[n]) {
                    Mindex.erase(Mindex.begin()+m);
                    Cindex.erase(Cindex.begin()+m);
                    Pindex.erase(Pindex.begin()+m);
                    Rindex.erase(Rindex.begin()+m);
                    Cyindex.erase(Cyindex.begin()+m);
                }
            }
        }
        for (m=0;m<cj_ref.size();m++) {
            for (n=0;n<aaa.Cindex.size();n++) {
                if (aaa.Cindex[n] == cj_ref[m]) {
                    aaa.Mindex.erase(aaa.Mindex.begin()+n);
                    aaa.Cindex.erase(aaa.Cindex.begin()+n);
                    aaa.Pindex.erase(aaa.Pindex.begin()+n);
                    aaa.Rindex.erase(aaa.Rindex.begin()+n);
                    aaa.Cyindex.erase(aaa.Cyindex.begin()+n);
                }
            }
        }
        t.clear();
        for (n=1;n<ci_ref.size();n++) {
            for (m=0;m<ci_ref.size();m++) {
                if (pi_ref[n] == ci_ref[m]) {
                    t.push_back(m);
                }
            }
        }
        f.clear();
        for (n=1;n<cj_ref.size();n++) {
            for (m=0;m<cj_ref.size();m++) {
                if (pj_ref[n] == cj_ref[m]) {
                    f.push_back(m);
                }
            }
        }
        pi_ref[0] = b;
        pj_ref[0] = a;
        ci_ref[0] = y;
        cj_ref[0] = x;

        for (n=1;n<ci_ref.size();n++) {
            ci_ref[n] = ci_ref[n-1] + 1;
            pi_ref[n] = ci_ref[t[n-1]];
        }
        for (n=1;n<cj_ref.size();n++) {
            cj_ref[n] = cj_ref[n-1] + 1;
            pj_ref[n] = cj_ref[f[n-1]];
        }

        t.clear();
        f.clear();
        for (n=1;n<Cindex.size();n++) {
            for (m=0;m<Cindex.size();m++) {
                if (Pindex[n] == Cindex[m]) {
                    t.push_back(m);
                }
            }
        }
        ref=-1;
        for (n=1;n<Cindex.size();n++) {
            if (Cindex[n] != Cindex[n-1]+1) {
                ref = n;
                Cindex[ref] = Cindex[ref-1] + cj_ref.size() + 1;
                break;
            }
        }
        if (ref!=-1) {
            for (n=ref+1;n<Cindex.size();n++) {
                Cindex[n] = Cindex[n-1] + 1;
            }

            for (n=1;n<Cindex.size();n++) {
                Pindex[n] = Cindex[t[n-1]];
            }
        }
        t.clear();
        f.clear();
        for (n=1;n<aaa.Cindex.size();n++) {
            for (m=0;m<aaa.Cindex.size();m++) {
                if (aaa.Pindex[n] == aaa.Cindex[m]) {
                    t.push_back(m);
                }
            }
        }
        ref=-1;
        for (n=1;n<aaa.Cindex.size();n++) {
            if (aaa.Cindex[n] != aaa.Cindex[n-1]+1) {
                ref = n;
                aaa.Cindex[ref] = aaa.Cindex[ref-1] + ci_ref.size() + 1;
                break;
            }
        }
        if (ref!=-1) {
            for (n=ref+1;n<aaa.Cindex.size();n++) {
                aaa.Cindex[n] = aaa.Cindex[n-1] + 1;
            }

            for (n=1;n<aaa.Cindex.size();n++) {
                aaa.Pindex[n] = aaa.Cindex[t[n-1]];
            }
        }
    }
    n=0;
    ref = 1;
    while(n<cj_ref.size()) {
        for (m=0;m<Cindex.size();m++) {
            if (Cindex[m] > cj_ref[n] && Cindex[m-1] < cj_ref[n]) {
                Cindex.insert(Cindex.begin() + m, cj_ref[n]);
                Pindex.insert(Pindex.begin() + m, pj_ref[n]);
                Mindex.insert(Mindex.begin() + m, mj_ref[n]);
                Rindex.insert(Rindex.begin() + m, rj_ref[n]);
                Cyindex.insert(Cyindex.begin()+m,cyj.at(n));
                n++;
                break;
            } else if (Cindex[Cindex.size()-1] < cj_ref[n]) {
                Cindex.push_back(cj_ref[n]);
                Pindex.push_back(pj_ref[n]);
                Mindex.push_back(mj_ref[n]);
                Rindex.push_back(rj_ref[n]);
                Cyindex.push_back(cyj.at(n));
                n++;
                break;
            }
        }
    }
    m=0;
    ref = 1;
    while (m<ci_ref.size()) {
        for (n=0;n<aaa.Cindex.size();n++) {
            if (aaa.Cindex[n] > ci_ref[m] && aaa.Cindex[n-1] < ci_ref[m]) {
                ref=n;
                aaa.Cindex.insert(aaa.Cindex.begin() + ref, ci_ref[m]);
                aaa.Pindex.insert(aaa.Pindex.begin() + ref, pi_ref[m]);
                aaa.Mindex.insert(aaa.Mindex.begin() + ref, mi_ref[m]);
                aaa.Rindex.insert(aaa.Rindex.begin() + ref, ri_ref[m]);
                aaa.Cyindex.insert(aaa.Cyindex.begin()+ref,cyi.at(m));
                m++;
                break;
            } else if (aaa.Cindex[aaa.Cindex.size()-1] < ci_ref[m]) {
                ref=n;
                aaa.Cindex.push_back(ci_ref[m]);
                aaa.Pindex.push_back(pi_ref[m]);
                aaa.Mindex.push_back(mi_ref[m]);
                aaa.Rindex.push_back(ri_ref[m]);
                aaa.Cyindex.push_back(cyi.at(m));
                m++;
                break;
            }
        }
    }
    reset();
	aaa.reset();
    vector<int>().swap(t);
    vector<int>().swap(f);
    vector<int>().swap(ci_ref);
    vector<int>().swap(mi_ref);
    vector<int>().swap(ri_ref);
    vector<int>().swap(cj_ref);
    vector<int>().swap(mj_ref);
    vector<int>().swap(rj_ref);
	vector<int>().swap(prj);
	vector<int>().swap(pri);
    return 1;
}


void MOLECULE::read(string a) {
	ifstream inf((a+".mds").c_str());
	string b;
	int i,j,k;
	while (!inf.eof()) {
		b="";
		inf>>b>>ws;
		if (b=="natom") inf>>j>>ws;
		else if (b=="Pindex") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				Pindex.push_back(k);
			}
		}
		else if (b=="Cindex") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				Cindex.push_back(k);
			}
		}
		else if (b=="Cyindex") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				Cyindex.push_back(k);
			}
		}
		else if (b=="Rindex") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				Rindex.push_back(k);
			}
		}
		else if (b=="if_circle") {
			inf>>k>>ws;
			if_circle=k;
		}
		else if (b=="Mindex") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				Mindex.push_back(k);
			}
		}
		else if (b=="Protect") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
			}
		}
	}
	return;
}

int MOLECULE::ring(int pt1,int pt2) {
    reset();
    int i,j,k,n=0,m=0,x,y;
    int b_pos[2][2];
    int b_circle;
	k=1;
	for (i=0;i<data.a.at(Mindex.at(pt1)).norder;i++) {
		if (Bindex[pt1][i]==1) k=0;
	}
	if (k) return 0;
	k=1;
	for (i=0;i<data.a.at(Mindex.at(pt2)).norder;i++) {
		if (Bindex[pt2][i]==1) k=0;
	}
	if (k) return 0;
    		i=pt1;
            j=pt2;
            if ((Pindex[i] != Pindex[j]) && (i!=j) && Pindex[i]!=Cindex[j] && Pindex[j]!=Cindex[i]&&Cyindex[i]==0&&Cyindex[j]==0) {
                for (x=0;x<data.a[Mindex[i]].norder;x++) {
                    if (Bindex[i][x]==1) {
                        for (y=0;y<data.a[Mindex[j]].norder;y++) {
                            if (Bindex[i][x] == Bindex[j][y]) {
                                b_pos[0][0]=x;
                                b_pos[0][1]=y;
                                b_pos[1][0]=i;
                                b_pos[1][1]=j;
                                b_circle=Bindex[i][x];
                                Bindex[i][x]=Bindex[j][y]=0;
                                Cyindex[i]=Cyindex[j]=if_circle+1;
                                m=1;
                                if_circle++;
                                break;
                            }
                        }
                        if (m==1) break;
                    }
                }
            }
			else {
				ofstream log("molecule.log",ios::app);
				log<<"Warning : this molecule "<<smiles<<" cannot be cyclized at point "<<pt1<<" and "<<pt2<<"."<<endl;
				log.close();
				return 0;
			}
	reset();
    return 1;
}

void MOLECULE::mds2smi()
{
	int j,q,n,w,i=Pindex.size(),t,tt,M;
	molesmi="";
    t=0;
    for (q=0;q<i;q++) {
        j=Mindex.at(q);
        j=data.a[j].nbond;
        t+=j;
    }
    if (if_circle) t+=2*if_circle;
    vector<char> x(1024,' ');

    if (if_circle) {
        n=0;
        vector<int> pos,m,num;
        char tmp[10];
        for (q=0;q<i;q++) {
            if (Cyindex.at(q) != 0) {
                m.push_back(q);
                pos.push_back(atomsmi[q][0]);
                num.push_back(Cyindex.at(q));
            }
        }
        for (q=0;q<pos.size();q++) {
            if (q==(pos.size()-1)) break;
            for (j=0;j<(pos.size()-q);j++) {
                if (j==(pos.size()-1)) break;
                if (pos[j]>pos[j+1]) {
                    swap(pos[j],pos[j+1]);
                    swap(m[j],m[j+1]);
                    swap(num[j],num[j+1]);
                }
            }
    }
        for (w=0;w<m.size();w++) {
            if (data.a[Mindex[m[w]]].chg) {
                if (num[w]>10) {
                    t=atomsmi[m[w]][5]+1;
                    tt=2;
                    sprintf(tmp,"%d%d",num[w]);
                    x[atomsmi[m[w]][5]+1]=tmp[0];
                    x[atomsmi[m[w]][5]+2]=tmp[1];
                }
                else {
                    tt=1;
                    t=atomsmi[m[w]][5]+1;
                    sprintf(tmp,"%d",num[w]);
                    x[atomsmi[m[w]][5]+1]=tmp[0];
                }
            }
            else {
                if (num[w]>10) {
                    t=atomsmi[m[w]][0]+1;
                    tt=2;
                    sprintf(tmp,"%d%d",num[w]);
                    x[atomsmi[m[w]][0]+1]=tmp[0];
                    x[atomsmi[m[w]][0]+2]=tmp[1];
                }
                else {
                    tt=1;
                    t=atomsmi[m[w]][0]+1;
                    sprintf(tmp,"%d",num[w]);
                    x[atomsmi[m[w]][0]+1] = tmp[0];
                }
            }
            for (q=0;q<i;q++) {
                j=Mindex.at(q);
                j=data.a[j].nbond;
                for (n=0;n<j;n++) {
                    if (atomsmi[q][n] >= t) atomsmi[q][n] = atomsmi[q][n]+tt;
                }
            }
        }
        pos.clear();
        m.clear();
        num.clear();
        vector<int>().swap(pos);
        vector<int>().swap(m);
        vector<int>().swap(num);
    }
    n=0;
    chg=0;
    for (q=0;q<Cindex.size();q++) chg+=data.a[Mindex.at(q)].chg;
    t=0;
   int tmp=0;
    for (q=0;q<i;q++) {
        j = Mindex.at(q);
        M = Mindex.at(q);
        j=data.a[j].nbond;
        t=t+j;
        tmp=0;
        if (data.a[M].chg) {
            for (tt=0;tt<data.a[M].norder;tt++) {
                tmp+=Bindex[q][tt];
            }
            if (tmp==1) data.a[M].name[3]='1';
            else if (tmp==2) data.a[M].name[3]='2';
            else if (tmp==3) data.a[M].name[3]='3';
            else if (tmp==4) data.a[M].name[3]='4';
            else if (tmp==5) data.a[M].name[3]='5';
        }
        tmp=data.a[Mindex[q]].index+1;
        for (w=0;w<j;w++) {
            n=atomsmi[q][w];
            if (((w-tmp)%3)==0 && w>data.a[Mindex[q]].index) {
                if (atomsmi[q][w-2]==(n-2) && atomsmi[q][w-1]==(n-1)) {
                    x[n]=' ';
                    x[n-1]=' ';
                    x[n-2]=' ';
                }
                else x[n]=data.a[Mindex[q]].name[w];
            }
            else x[n]=data.a[Mindex[q]].name[w];
        }
    }
    if (if_circle) t=t+2*if_circle;
    w=0;
    for (j=0;j<x.size();j++) {
        if (x[j]!=' ') {
            molesmi.push_back(x[j]);
        }
    }
    return;
}

void MOLECULE::clear() {

	int n;
	molesmi.clear();
	for (n=0;n<1024;n++) {
		Bindex[n].clear();
		atomsmi[n].clear();
	}
	return;
}


void MOLECULE::reset() {
	int i;
	int j;
	int n,m;
	int b;
	int c;
	int x;
	int u;
	int M;
	int smindex,t[2],s[2];
	clear();
	for (i=0;i<Cindex.size();i++) {
		b = Mindex.at(i);
		c = Pindex.at(i);
		x = Rindex.at(i);
		for (j=0;j<data.a[b].norder;j++) {
			Bindex[i].push_back(data.a[b].order[j]);
		}
		u=0;
		if (c==0) {
			for (j=0;j<data.a[b].nbond;j++) {
				atomsmi[i].push_back(j);
			}
		}
		else if (c>0) {
			M = Mindex.at(c-1);
			for (j=0;j<data.a[M].norder;j++) {
				if (Bindex[c-1].at(j) == x) {
					for (n=0;n<data.a[b].norder;n++) {
						if (Bindex[i].at(n) == x) {
							Bindex[i].at(n) = 0;
							Bindex[c-1].at(j) = 0;
							u=1;
							break;
						}
					}
					if (u==1) break;
				}
			}
			smindex=0;
			if (u) smindex=data.a[M].index+j*3;
			for (j=0;j<data.a[b].nbond;j++) {
				atomsmi[i].push_back(j+atomsmi[c-1][smindex]+1);
			}
			for (j=0;j<i;j++) {
				for (n=0;n<data.a[Mindex.at(j)].nbond;n++) {
					if (atomsmi[j].at(n) > atomsmi[c-1].at(smindex)) {
						atomsmi[j].at(n) = atomsmi[j].at(n)+data.a[b].nbond;
					}  
				}
			}
		}
	}
	m=0;
    for (i=0;i<Cyindex.size();i++) {
        if (Cyindex.at(i)>m) {
            if (Cyindex.at(i)>10) {
                t[0]=Cyindex.at(i)/10;
                t[1]=Cyindex.at(i)%10;
                if (t[0]>m && t[0]>t[1]) m=t[0];
                else if (t[1]>t[0] && t[1]>m) m=t[1];
            }
            else m=Cyindex.at(i);
        }
    }
	if (if_circle) {
		vector<int> cyclic;
		cyclic.resize(m);
		for (i=0;i<m;i++) cyclic[i]=0;
		for (i=0;i<Mindex.size();i++) {
            if (Cyindex.at(i))  {
                if (Cyindex.at(i)>10) {
                    j=Cyindex.at(i)/10;
                    n=Cyindex.at(i)%10;
                    cyclic.at(j-1)+=1;
                    cyclic.at(n-1)+=1;
                }
                else cyclic.at(Cyindex.at(i)-1)+=1;
            }
        }
		m=-1;
		for (i=0;i<cyclic.size();i++) {
            if (cyclic.at(i)==2) {
                u=0;
                for (j=0;j<Mindex.size();j++) {
                    for (n=j+1;n<Mindex.size();n++) {
                        if (Cyindex.at(j)>10) {
                            t[0]=Cyindex.at(j)/10;
                            t[1]=Cyindex.at(j)%10;
                        }
                        else {
                            t[0]=t[1]=Cyindex.at(j);
                        }
                        if (Cyindex.at(n)>10) {
                            s[0]=Cyindex.at(n)/10;
                            s[1]=Cyindex.at(n)%10;
                        }
                        else {
                            s[0]=s[1]=Cyindex.at(n);
                        }
                        if ( (s[0]==t[0] || s[0]==t[1] || s[1]==t[0] || s[1]==t[1]) && (s[0]==(i+1) || s[1]==(i+1))) {
                            for (m=0;m<data.a[Mindex.at(j)].norder;m++) {
                                for (x=0;x<data.a[Mindex.at(n)].norder;x++) {
                                    if (Bindex[j].at(m)==1 && Bindex[n].at(x)==1) {
                                        Bindex[j].at(m)=0;
                                        Bindex[n].at(x)=0;
                                        u=1;
                                        break;
                                    }
                                }
                                if (u) break;
                            }
                            if (u==0) {
                                u=1;
                                Cyindex.at(j)=0;
                                Cyindex.at(n)=0;
}
                            if (u) break;
                        }
                        if (u) break;
                    }
                    if (u) break;
                }
            } else {
                for (j=0;j<Mindex.size();j++) {
                    if (Cyindex.at(j)==(i+1)) {
                        Cyindex.at(j)=0;
                        break;
                    }
                }
            }
        }
    }
    n=0;
    if (if_circle) {
        m=0;
        for (i=0;i<Cyindex.size();i++) {
            if (Cyindex.at(i)>m) {
                if (Cyindex.at(i)>10) {
                    j=Cyindex.at(i)/10;
                    n=Cyindex.at(i)%10;
                    if (j>n && j>m) m=j;
                    else if (n>j&&n>m) m=n;
                }
                else m=Cyindex.at(i);
            }
        }
        if_circle=m;
    } 
	return;
}

void MOLECULE::init() {
	int i=0,j=0,k=0;
	string tmp;
	double x,y,z;
	smi2gjf();
	if (atm!=NULL) delete [] atm;
	string temp="";
	ifstream inf((smiles+".gjf").c_str());
    getline(inf,temp);
    inf>>ws;
    getline(inf,temp);
	inf>>ws;
	while (!inf.eof()) {
		inf>>tmp>>ws>>x>>ws>>y>>ws>>z>>ws;
		i++;
	}
	inf.close();
	natom=i;
	inf.open((smiles+".gjf").c_str());
	k=0;
	atm = new DEATOM [natom];
    getline(inf,temp);
    inf>>ws;
    getline(inf,temp);
    inf>>ws;
	while (!inf.eof()) {
			inf>>atm[k].name>>ws>>atm[k].x[0]>>ws>>atm[k].x[1]>>ws>>atm[k].x[2]>>ws;
			atm[k].find_r();
			k++;
			if (k==natom) break;
	}
	dist = new double *[natom];
	connect = new int *[natom];
	order= new int *[natom];
	for (i=0;i<natom;i++) {
		order[i]=new int [natom];
		connect[i]=new int [natom];
		dist[i]=new double [natom];
	}
	for (i=0;i<natom;i++) {
		dist[i][i]=0.0;
	}
	k=0;
	for (i=0;i<natom;i++) {
		if (atm[i].name != "H") {
			k++;
		}
	}
	for (i=0;i<k;i++) {
		Cindex.push_back(i+1);
		Cyindex.push_back(0);
	}
	inf.close();
	if_circle=0;
	cal_r();
	check_bnd();
	return;
}

void MOLECULE::smi2gjf() {
    system("pwd > pwddir.txt");
    ifstream PWD("pwddir.txt");
    string pwd="";
    PWD >> pwd >> ws;
    PWD.close();

	ofstream out((pwd+"/"+smiles+".smi").c_str());
	out << smiles << endl;
	out.close();
	PWD.open((smiles+".gjf").c_str());
	if (!PWD.is_open()){
		ofstream out("molecule.log",ios::app);
		out<<"Error : you don't include the gjf file of molecule, "<<smiles<<endl;
		out.close();
		exit(1);
	}
	return;
}

void MOLECULE::cal_r() {
	int i,j,k;
	double r;
	for (i=0;i<natom;i++) {
		for (j=i+1;j<natom;j++) {
			r = (atm[i].x[0]-atm[j].x[0])*(atm[i].x[0]-atm[j].x[0])+(atm[i].x[1]-atm[j].x[1])*(atm[i].x[1]-atm[j].x[1])+(atm[i].x[2]-atm[j].x[2])*(atm[i].x[2]-atm[j].x[2]);
			r = sqrt(r);
			dist[i][j] = dist[j][i] = r;
		}
	}
	return;
}
void MOLECULE::check_bnd(int aaa) {
    int i,j,k,sum[natom];
    double r,tol=1.15;
    for (i=0;i<natom;i++) {
        connect[i][i]=0;
    }
    for (i=0;i<natom;i++) {
        for (j=i+1;j<natom;j++) {
            r=tol*(atm[i].r_bnd+atm[j].r_bnd);
            if (dist[i][j]<=r) connect[i][j]=connect[j][i]=1;
            else connect[i][j]=connect[j][i]=0;
        }
    }
	for (i=0;i<natom;i++) {
        sum[i]=0;
        for (j=0;j<natom;j++) sum[i]+=connect[i][j];
    }

    for (i=0;i<natom;i++) {
        for (j=0;j<natom;j++) order[i][j]=connect[i][j];
    }
    k=0;
    for (i=0;i<natom;i++) {
		k=0;
        while (1) {
            if (atm[i].nbnd > sum[i]) {
                for (j=i+1;j<natom;j++) {
                    if ( atm[i].nbnd != sum[i] && connect[i][j] !=0 && atm[j].nbnd!=sum[j]) {
                        if ((atm[j].nbnd-sum[j])>0) {
                            order[i][j]+=1;
                            sum[i]+=1;
                            sum[j]+=1;
                            order[j][i] = order[i][j];
                        }
                    }
                }
				k++;
            }
            if (atm[i].nbnd==sum[i]) break;
			else if (k>100) break;
        }
    }

    for (i=0;i<natom;i++) {
        for (j=0;j<natom;j++) {
            connect[i][j]=order[i][j];
        }
    }
    return;
}

void MOLECULE::check_bnd() {
	int i,j,k,sum[natom],chgg;
	string tmp;
	double r,tol=1.15;
	ifstream inf((smiles+".gjf").c_str());
	getline(inf,tmp);
	inf>>ws;
	inf>>chgg>>ws>>k>>ws;
	inf.close();
	for (i=0;i<natom;i++) {
		connect[i][i]=0;
	}
	for (i=0;i<natom;i++) {
		for (j=i+1;j<natom;j++) {
			r=tol*(atm[i].r_bnd+atm[j].r_bnd);
			if (dist[i][j]<=r) connect[i][j]=connect[j][i]=1;
			else connect[i][j]=connect[j][i]=0;
		}
	}
    for (i=0;i<natom;i++) {
        if (atm[i].name!="C"&&atm[i].name!="H"&&chgg) {
            k=0;
            k=rd_nps(i);
			if (k==-1) {}
            else if (k!=atm[i].nbnd) atm[i].nbnd=k;
        }
    }

	for (i=0;i<natom;i++) {
		sum[i]=0;
		for (j=0;j<natom;j++) sum[i]+=connect[i][j];
	}

	for (i=0;i<natom;i++) {
		for (j=0;j<natom;j++) order[i][j]=connect[i][j];
	}
	k=0;
    for (i=0;i<natom;i++) {
        while (1) {
			k=0;
            if (atm[i].nbnd > sum[i]) {
                for (j=i+1;j<natom;j++) {
                    if ( atm[i].nbnd != sum[i] && connect[i][j] !=0 && atm[j].nbnd!=sum[j]) {
                        if ((atm[j].nbnd-sum[j])>0) {
                            order[i][j]+=1;
                            sum[i]+=1;
                            sum[j]+=1;
                            order[j][i] = order[i][j];
                        }
                    }
                }
				k++;
            }
            if (atm[i].nbnd==sum[i]) break;
			else if (k>100) break;
        }
    }
	for (i=0;i<natom;i++) {
		for (j=0;j<natom;j++) {
			connect[i][j]=order[i][j];
		}
	}
	return;
}
void MOLECULE::smi2cod() {
	int i,j,n,t,m,u;
	vector<int> k;
	for (i=0;i<natom;i++) {
		k.clear();
		if (atm[i].name == "H") continue; 
		else if (atm[i].name == "C") { 
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size() == 4) Mindex.push_back(1);
			else if (k.size() == 3) {
				Mindex.push_back(2);
			}
			else if (k.size()==2) {
				for (n=0;n<k.size();n++) {
					if (k[n]==3 || k[n]==1) {
						Mindex.push_back(3);
						break;
					} else if (k[n]==2) {
						Mindex.push_back(4);
						break;
					}
				}
			}
		}
		else if (atm[i].name=="O") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==1) {
				if (k[0]==2) Mindex.push_back(6);
				else if (k[0]==1) Mindex.push_back(29);
			}
			else if (k.size()==2) {
				for (j=0;j<natom;j++) {
					if (connect[i][j] != 0) {
						if (atm[j].name=="H") {
							Mindex.push_back(5);
							n=0;
							break;
						} else n=1;
					}
				}
				if (n==1) Mindex.push_back(5);
			}
		}
		else if (atm[i].name=="N") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==4) Mindex.push_back(15);
			if (k.size()==3) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==4) Mindex.push_back(16);
				else if (t==3) Mindex.push_back(7);
			}
			if (k.size()==2) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==3) Mindex.push_back(8);
				else if (t==2) Mindex.push_back(35);
			}
			if (k.size()==1) Mindex.push_back(9);
		}
		else if (atm[i].name=="P") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==6) Mindex.push_back(30);
			if (k.size()==5) Mindex.push_back(31);
			if (k.size()==4) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==4) Mindex.push_back(17);
				else if (t==5) Mindex.push_back(32);
			}
			if (k.size()==3) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==4) Mindex.push_back(18);
				else if (t==3) Mindex.push_back(21);
			}
			if (k.size()==2) Mindex.push_back(22);
			if (k.size()==1) Mindex.push_back(23);
		}
		else if (atm[i].name=="F") {
            for (j=0;j<natom;j++) {
                if (connect[i][j] != 0) k.push_back(connect[i][j]);
            }
            if (k.size()) Mindex.push_back(11);
            else Mindex.push_back(24);
        }
        else if (atm[i].name=="Cl") {
            for (j=0;j<natom;j++) {
                if (connect[i][j] != 0) k.push_back(connect[i][j]);
            }
            if (k.size()) Mindex.push_back(12);
            else Mindex.push_back(25);
        }
        else if (atm[i].name=="Br") {
            for (j=0;j<natom;j++) {
                if (connect[i][j] != 0) k.push_back(connect[i][j]);
            }
            if (k.size()) Mindex.push_back(13);
            else Mindex.push_back(26);
        }
        else if (atm[i].name=="I") {
            for (j=0;j<natom;j++) {
                if (connect[i][j] != 0) k.push_back(connect[i][j]);
            }
            if (k.size()) Mindex.push_back(14);
            else Mindex.push_back(27);
        }

		else if (atm[i].name=="S") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==2) Mindex.push_back(19);
			else if (k.size()==1) Mindex.push_back(20);
			else if (k.size()==3) Mindex.push_back(34);
			else if (k.size()==4) Mindex.push_back(33);
		} 
	}
	Rindex.push_back(-1);
	Pindex.push_back(0);
	vector<int> c_ref;
	for (i=0;i<natom;i++) {
		if (atm[i].name != "H") {
			c_ref.push_back(i);
		}
	}
	for (i=0;i<natom;i++) {
		if (atm[i].name=="H") {
			for (j=0;j<natom;j++) {
				order[i][j]=0;
				order[j][i]=0;
			}
		}
	}
	int tmp[Cindex.size()];
	n=0;
	for (i=0;i<natom;i++) {
		if (atm[i].name!="H") {
			tmp[n]=i;
			n++;
		}
	}

	for (i=1;i<natom;i++) {
		t=0;
		m=0;
		if (atm[i].name!="H") {
			for (j=i;j>=0;j--) {
				if (order[i][j]>m) {
					m=order[i][j];
				}
			}
			for (j=i;j>=0;j--) {
				if (order[i][j]==m) {
					for (n=0;n<Cindex.size();n++) {
						if (tmp[n]==j) {
							Rindex.push_back(order[i][j]);
							Pindex.push_back(Cindex[n]);
							order[j][i]=order[i][j]=0;
							t=1;
							break;
						}
					}
				}
				if (t) break;
			}
		}
	}
	n=1;
	for (i=0;i<natom;i++) {
		for (j=i+1;j<natom;j++) {
			if (atm[j].name!="H" && order[i][j] != 0) {
				Cyindex.at(i)=Cyindex.at(i)*10+n;
				Cyindex.at(j)=Cyindex.at(j)*10+n;
				order[i][j]=order[j][i]=0;
				if_circle=n;
				n++;
			}
		}
	}
	return;
}

void MOLECULE::report() {
	ofstream outf((smiles+".mds").c_str());
	int i;
	outf<<"natom "<<Cindex.size();
	outf<<endl;
	outf<<"Mindex ";
	for (i=0;i<Cindex.size();i++) outf<<Mindex[i]<<" ";
	outf<<endl;
	outf<<"Pindex ";
	for (i=0;i<Cindex.size();i++) outf<<Pindex[i]<<" ";
	outf<<endl;
	outf<<"Cindex ";
	for (i=0;i<Cindex.size();i++) outf<<Cindex[i]<<" ";
	outf<<endl;
	outf<<"Rindex ";
	for (i=0;i<Cindex.size();i++) outf<<Rindex[i]<<" ";
	outf<<endl;
	outf<<"Cyindex ";
	for (i=0;i<Cindex.size();i++) {
		outf<<Cyindex[i]<<" ";
	}
	outf<<endl;
	outf<<"if_circle "<<if_circle<<endl;
	Cindex.clear();
	Mindex.clear();
	Pindex.clear();
	Cyindex.clear();
	Rindex.clear();
	if_circle=0;
	outf.close();
	clean();
	return;
}

void MOLECULE::clean() {
	int i;
	if(dist==NULL) return;
	for (i=0;i<natom;i++) {
		delete[] dist[i];
		delete[] connect[i];
		delete[] order[i];
	}
	delete[] dist;
	delete[] connect;
	delete[] order;
	delete[] atm;
	dist=NULL;
	connect=NULL;
	order=NULL;
	atm=NULL;
	return;
}


void MOLECULE::outmds(ostream &outs) {
	int i,j,k;
    outs<<"natom "<<Cindex.size();
    outs<<endl;
    outs<<"Mindex ";
    for (i=0;i<Cindex.size();i++) outs<<Mindex[i]<<" ";
    outs<<endl;
    outs<<"Pindex ";
    for (i=0;i<Cindex.size();i++) outs<<Pindex[i]<<" ";
    outs<<endl;
    outs<<"Cindex ";
    for (i=0;i<Cindex.size();i++) outs<<Cindex[i]<<" ";
    outs<<endl;
    outs<<"Rindex ";
    for (i=0;i<Cindex.size();i++) outs<<Rindex[i]<<" ";
    outs<<endl;
    outs<<"Cyindex ";
    for (i=0;i<Cindex.size();i++) {
        outs<<Cyindex[i]<<" ";
    }
    outs<<endl;
    outs<<"if_circle "<<if_circle<<endl<<endl;
    return;
}

int MOLECULE::combine(MOLECULE &B, int k,int p) {
    int i,j,tmp,counter,x;
    x=1;
    for (i=0;i<data.a[Mindex[k]].norder;i++) {
        for (j=0;j<B.data.a[B.Mindex[p]].norder;j++) {
            if (Bindex[k].at(i)>=1 && Bindex[k].at(i)==B.Bindex[p].at(j)) x=0;
        }
    }
    if (x) {
        reset();
        mds2smi();
		ofstream log("molecule.log",ios::app);
		log<<"Warning : "<<smiles<<" cannot combine with "<<B.smiles<<" at point "<<k<<" and "<<p<<" because of different bond order."<<endl;
        return 0;
    }
    if (p==0) ;
    else {
        B.recode(p);
		B.reset();
    }
    if (if_circle && B.if_circle) {
        for (i=0;i<B.Cindex.size();i++) {
            if (B.Cyindex.at(i)>0) B.Cyindex.at(i)+=if_circle;
        }
    }
    counter=0;
    while (1) {
        tmp=0;
        for (i=0;i<data.a[Mindex.at(k)].norder;i++) {
            if (Bindex[k].at(i)>0) {
                j=i; // the position
                tmp=Bindex[k].at(i); // the order
                break;
            }
        }
		x=1;
		for (i=0;i<B.data.a.at(B.Mindex.at(0)).norder;i++) {
			if (B.Bindex[0][i]==tmp) {
				x=0;
				break;
			}
		}
		if (x) {
			B.clear();
			B.empty();
			B.input();
			return 0;
		}
        if (tmp!=0) {
            x=1;
            for (i=0;i<B.data.a[B.Mindex[0]].norder;i++) {
                if (B.Bindex[0][i]==tmp) {
                    x=0;
                    break;
                }
            }
            if (1) {
                B.Rindex[0]=tmp;
				int tmpp=Cindex.size();
                for (i=0;i<B.Cindex.size();i++) {
                    B.Cindex[i]+=tmpp;
                    if (i==0) B.Pindex[i]=Cindex[k];
                    else B.Pindex[i]+=tmpp;
                }
                for (i=0;i<B.Cindex.size();i++) {
                    Pindex.push_back(B.Pindex.at(i));
                    Cindex.push_back(B.Cindex.at(i));
                    Rindex.push_back(B.Rindex.at(i));
                    Cyindex.push_back(B.Cyindex.at(i));
                    Mindex.push_back(B.Mindex.at(i));
                }
                B.clear();
                B.empty();
                B.input();
                if_circle=B.if_circle+if_circle;
                reset();
                mds2smi();
                return 1;
            }

        }
        else if (tmp==0) {
            if (p) {
                B.clear();
                B.empty();
                B.input();
            }
            reset();
            mds2smi();
            return 0;
        }
    }
    return 1;
}


void MOLECULE::recode(int k) {
	if (smiles=="") smiles=molesmi;
	if (1) {
	empty();
	clear();
	init(k);
	check_bnd(k);
	smi2cod();
	report();
	read(smiles);
	reset();
	mds2smi();
	}
	return;
}
void MOLECULE::wipe() {
	vector<int>().swap(Cindex);
	vector<int>().swap(Pindex);
	vector<int>().swap(Cyindex);
	vector<int>().swap(Rindex);
	vector<int>().swap(Mindex);
	int i;
	for (i=0;i<1024;i++) {
		vector<int>().swap(Bindex[i]);
		vector<int>().swap(atomsmi[i]);
	}
	delete [] atm;
	molesmi.clear();
	smiles.clear();
	return;
}
void MOLECULE::init(int p) {
	int i=0,j=0,k=0,n=0,s=0,q=0;
	vector<string> n_tmp;
	string trash;
	MOLECULE tmp;
	double **d_tmp=NULL;
	double x,y,z;
	smi2gjf();
	string* temp=new string("a");
	tmp.smiles=smiles;
	tmp.init();
	int point[tmp.Cindex.size()];
	ifstream inf((smiles+".gjf").c_str());
	getline(inf,*temp);
	inf>>ws;
	getline(inf,*temp);
	inf>>ws;
	i=0;
	while (!inf.eof()) {
		inf>>trash>>x>>y>>z>>ws;
		i++;
	}
	inf.close();
	natom=i;
	inf.open((smiles+".gjf").c_str());
	atm=new DEATOM [natom];
	d_tmp=new double *[natom];
	for (i=0;i<natom;i++) {
		d_tmp[i]=new double [3];
	}
    getline(inf,*temp);
    inf>>ws;
    getline(inf,*temp);
    inf>>ws;
	k=0;
	while (!inf.eof()) {
			inf>>trash>>d_tmp[k][0]>>d_tmp[k][1]>>d_tmp[k][2];
			n_tmp.push_back(trash);
			k++;
			if (k==natom) break;
	}
	inf.close();
   	delete temp;
	k=0;
	for (i=0;i<natom;i++) {
		if (n_tmp[i]!="H") {
			point[k]=i;
			k++;
		}
	}

	n=point[p];
	atm[0].name=n_tmp[point[p]];
	for (j=0;j<3;j++) atm[0].x[j]=d_tmp[point[p]][j];
	if (atm[0].name!="C"&&atm[0].name!="H") atm[0].nbnd=tmp.atm[point[p]].nbnd;
	atm[0].find_r();
    point[p]=-1;
	q=1;
	j=0;
	s=0;
	while (1) {
		if (q==1) {
			j=0;
			for (i=0;i<natom;i++) {
				if (tmp.connect[i][n]==3) {j=i; break;}
				else if (tmp.connect[i][n]==2) {j=i; break;}
			}
			if (j>n) s=1;
			else s=0;
			if (s) {
				for (i=0;i<k;i++) {
					if (point[i]==j) {
						s=point[i];
						point[i]=-1;
						break;
					} else s=-2; 
				}
				if (s!=-2) {
					for (i=0;i<natom;i++) {
						//tmp.connect[i][n]=0;
						tmp.connect[n][i]=0;
					}
					atm[q].name=n_tmp[s];
					for (i=0;i<3;i++) atm[q].x[i]=d_tmp[s][i];
					atm[q].find_r();
					if (atm[q].name!="C"&&atm[q].name!="H") atm[q].nbnd=tmp.atm[s].nbnd;
					q++;
				}
			}
		}	

		for (i=n;i>=0;i--) {
			if (tmp.connect[i][n] != 0 && n_tmp[i] != "H") {
				for (j=0;j<k;j++) {
					if (point[j]==i) {
						s=point[j];
						point[j]=-1;
						break;
					} else s=-2;
				}
				if (s==-2) break;
				for (j=0;j<natom;j++) {
					tmp.connect[j][n]=0;
					tmp.connect[n][j]=0;
				}
				atm[q].name=n_tmp[s];
				for (j=0;j<3;j++) atm[q].x[j]=d_tmp[s][j];
				atm[q].find_r();
				if (atm[q].name!="H"&&atm[q].name!="C") atm[q].nbnd=tmp.atm[s].nbnd;
				q++;
				n=i;
				break;
			}
		}
		s++;
		if (s==natom) break;

	}
	for (i=0;i<k;i++) {
		if (point[i]!=-1) {
			atm[q].name=n_tmp[point[i]];
			for (j=0;j<3;j++) atm[q].x[j]=d_tmp[point[i]][j];
			atm[q].find_r();
			if (atm[q].name=="N"||atm[q].name=="S"||atm[q].name=="P") {atm[q].nbnd=tmp.atm[point[i]].nbnd;}
			q++;
		}
	}
	for (i=0;i<natom;i++) {
		delete [] tmp.dist[i];
		delete [] tmp.connect[i];
		delete [] tmp.order[i];
	}
	delete [] tmp.dist;
	delete [] tmp.connect;
	delete [] tmp.order;
	tmp.wipe();
	for (i=0;i<natom;i++) {
		if (q>=natom) break;
		if (n_tmp[i]=="H") {
			atm[q].name=n_tmp[i];
			for (j=0;j<3;j++) {
				atm[q].x[j]=d_tmp[i][j];
			}       
			atm[q].find_r();
			q++;
		}
	}

	for (i=0;i<natom;i++) {
		delete [] d_tmp[i];
	}
	delete [] d_tmp;
	d_tmp=NULL;

	dist = new double *[natom];
	connect = new int *[natom];
	order= new int *[natom];
	for (i=0;i<natom;i++) {
		order[i]=new int [natom];
		connect[i]=new int [natom];
		dist[i]=new double [natom];
	}
	for (i=0;i<natom;i++) {
		dist[i][i]=0.0;
	}
	k=0;
	for (i=0;i<natom;i++) {
		if (atm[i].name != "H") {
			k++;
		}
	}
	for (i=0;i<k;i++) {
		Cindex.push_back(i+1);
		Cyindex.push_back(0);
	}
	inf.close();
	if_circle=0;
	cal_r();
	return;
}

void MOLECULE::input() {
	empty();
	if (atm!=NULL) delete [] atm;
	init();
	smi2cod();
	report();
	read((smiles).c_str());
	reset();
	mds2smi();
	return;
}

void MOLECULE::empty() {
	Cindex.clear();
	Pindex.clear();
	Mindex.clear();
	Rindex.clear();
	Cyindex.clear();
	molesmi="";
	for (int i=0;i<1024;i++) {
		Bindex[i].clear();
		atomsmi[i].clear();
	}
	return;
}

int MOLECULE::rd_nps(int a) {
    system("pwd > pwddir.txt");
    ifstream PWD("pwddir.txt");
    string pwd="";
    PWD >> pwd >> ws;
    PWD.close();
	ifstream inf((smiles+".mol").c_str());
	if (!inf.is_open()) {
		ofstream out("molecule.log",ios::app);
		out<<"Warning : you don't include the mol file of this molecule, "<<smiles<<". It might cause the error of ionic molecule."<<endl;
		out.close();
		return -1;
	}
	string tmp;
	int i,j,k,l=0,s=0,m,n,p,q;
	a++;
	inf>>ws;
	getline(inf,tmp);
	inf>>ws;
	getline(inf,tmp);
	for (i=0;i<natom;i++) getline(inf,tmp);
	for (i=0;i<natom;i++) {
		inf>>j>>k>>l>>m>>n>>p>>q>>ws;
		if (j==a || k==a) s+=l;
	}
	inf.close();
	return s;
}
int MOLECULE::add(int pt, int id) {
	int i,j,n=0;
	for (i=0;i<data.a[Mindex.at(pt)].norder;i++) {
		for (j=0;j<data.a[id].norder;j++) {
			if (Bindex[pt][i]==data.a[id].order[j]) {
                Pindex.push_back(Cindex.at(pt));
                Cindex.push_back(Cindex.size()+1);
                Mindex.push_back(id);
                Rindex.push_back(Bindex[pt][i]);
                Cyindex.push_back(0);
				n=1;
				break;
			}
		}
		if (n) break;
	}
	if (n==0) {
		ofstream log("molecule.log",ios::app);
		log<<"Warning : "<<"for "<<smiles<<", it can't excute addition at this point "<<pt<<" with this atom (id), "<<id<<"."<<endl;
		log.close();
		reset();
		return 0;
	}
	reset();
	return 1;
}

void MOLECULE::subtract(int n) {
	int i,j,k;
	vector<int> ref;
	ref.push_back(Cindex.at(n));
	if (n==0) {
		ofstream log("molecule.log",ios::app);
		log<<"Warning : This molecule would be empty because of point, "<<n<<endl;
		log.close();
	}
	for (i=n;i<Cindex.size();i++) {
		if (Cindex[i]>n) {
			j=ref.size();
			for (k=0;k<j;k++) {
				if (Pindex[i]==ref[k]) {
					ref.push_back(Cindex.at(i));
				}
			}
		}
	}
	for (i=0;i<ref.size();i++) {
		for (j=0;j<Cindex.size();j++) {
			if (ref[i]==Cindex.at(j)) {
				Mindex.erase(Mindex.begin()+j);
				Cindex.erase(Cindex.begin()+j);
				Pindex.erase(Pindex.begin()+j);
				Rindex.erase(Rindex.begin()+j);
				Cyindex.erase(Cyindex.begin()+j);	
			}
		}
	}
	for (i=0;i<Cindex.size();i++) {
		if (Cindex.at(i)!=(i+1)) {
			j=Cindex.at(i);
			Cindex.at(i)=i+1;
			for (k=i;k<Cindex.size();k++) {
				if (Pindex.at(k)==j) Pindex.at(k)=Cindex.at(i);
			}
		}
	}
	reset();
	return;
}
int MOLECULE::exchange(int n,int id,int id2,int bond) {
    int i,j,k,M,P,C,R,bd[3];
	reset();
    vector<int> tmp;
    for (i=0;i<data.a[Mindex.at(n)].norder;i++) tmp.push_back(Bindex[n].at(i));
    M=Mindex.at(n);
    P=Pindex.at(n);
    C=Cindex.at(n);
    R=Rindex.at(n);
		if (R+bond<1 || R+bond>3) {
			ofstream log("molecule.log",ios::app);
			log<<"Warning : for "<<smiles<<", it can't be exchanged at point, "<<n<<" to become as "<<id<<" and "<<id2<<" with add "<<bond<<" order."<<endl;
			log.close();
			reset();
			return 0;
		}
        if (R==2) {
            bd[0]=bd[1]=bd[2]=0;
            for (i=0;i<tmp.size();i++) {
                if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
            }
            tmp.clear();
            for (i=0;i<4;i++) bd[i]=data.a[M].bd[i]-bd[i];
			if (bond<0) { bd[0]+=1;bd[1]-=1;}
			else if (bond>0) {bd[2]+=1;bd[1]-=1;}
            k=0;j=0;i=id;
            if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[M].chg==data.a[i].chg) {Mindex.at(n)=i;Rindex.at(n)=R+bond;}
			else {
				reset();
            	ofstream log("molecule.log",ios::app);
            	log<<"Warning : for "<<smiles<<", it can't be exchanged at point, "<<n<<" to become as "<<id<<" and "<<id2<<" with add "<<bond<<" order."<<endl;
            	log.close();
				return 0;
			}
            bd[0]=bd[1]=bd[2]=0;
            for (i=0;i<data.a[Mindex[P-1]].norder;i++) tmp.push_back(Bindex[P-1].at(i));
            for (i=0;i<tmp.size();i++) {
                if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
            }
            tmp.clear();
            for (i=0;i<4;i++) bd[i]=data.a[Mindex.at(P-1)].bd[i]-bd[i];
            if (bond<0) {bd[0]+=1;bd[1]-=1;}
			else if (bond>0) {bd[2]+=1;bd[1]-=1;}
            k=0;j=0;i=id2;
            if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[Mindex[P-1]].chg==data.a[i].chg) {Mindex.at(P-1)=i;Rindex.at(n)=R+bond;}
			else {
            	ofstream log("molecule.log",ios::app);
            	log<<"Warning : for "<<smiles<<", it can't be exchanged at point, "<<n<<" to become as "<<id<<" and "<<id2<<" with add "<<bond<<" order."<<endl;
            	log.close();
				Rindex.at(n)=R;
				Mindex.at(n)=M;
				reset();
				return 0;
			}
		}
		else if (R==3) {
            bd[0]=bd[1]=bd[2]=0;
            for (i=0;i<tmp.size();i++) {
                if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
            }
            tmp.clear();
            for (i=0;i<4;i++) bd[i]=data.a[M].bd[i]-bd[i];
            if (bond==-1) { bd[1]+=1;bd[2]-=1;}
            else if (bond==-2) {bd[0]+=1;bd[2]-=1;}
            k=0;j=0;i=id;
            if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[M].chg==data.a[i].chg) {Mindex.at(n)=i;Rindex.at(n)=R+bond;}
            else {
            	ofstream log("molecule.log",ios::app);
            	log<<"Warning : for "<<smiles<<", it can't be exchanged at point, "<<n<<" to become as "<<id<<" and "<<id2<<" with add "<<bond<<" order."<<endl;
            	log.close();
                reset();
                return 0;
            }
            bd[0]=bd[1]=bd[2]=0;
            for (i=0;i<data.a[Mindex[P-1]].norder;i++) tmp.push_back(Bindex[P-1].at(i));
            for (i=0;i<tmp.size();i++) {
                if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
            }
            tmp.clear();
            for (i=0;i<4;i++) bd[i]=data.a[Mindex.at(P-1)].bd[i]-bd[i];
            if (bond==-1) {bd[1]+=1;bd[2]-=1;}
            else if (bond==-2) {bd[0]+=1;bd[2]-=1;}
            k=0;j=0;i=id2;
            if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[Mindex[P-1]].chg==data.a[i].chg) {Mindex.at(P-1)=i;Rindex.at(n)=R+bond;}
            else {
            	ofstream log("molecule.log",ios::app);
            	log<<"Warning : for "<<smiles<<", it can't be exchanged at point, "<<n<<" to become as "<<id<<" and "<<id2<<" with add "<<bond<<" order."<<endl;
            	log.close();
				Rindex.at(n)=R;
                Mindex.at(n)=M;
                reset();
                return 0;
            }
        }
        else if (R==1) {
            bd[0]=bd[1]=bd[2]=0;
            for (i=0;i<tmp.size();i++) {
                if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
            }
            tmp.clear();
            for (i=0;i<4;i++) bd[i]=data.a[M].bd[i]-bd[i];
            if (bond==1) { bd[1]+=1;bd[0]-=1;}
            else if (bond==2) {bd[2]+=1;bd[0]-=1;}
            k=0;j=0;i=id;
            if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[M].chg==data.a[i].chg) {Mindex.at(n)=i;Rindex.at(n)=R+bond;}
            else {
                reset();
                return 0;
            }
            bd[0]=bd[1]=bd[2]=0;
            for (i=0;i<data.a[Mindex[P-1]].norder;i++) tmp.push_back(Bindex[P-1].at(i));
            for (i=0;i<tmp.size();i++) {
                if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
            }
            tmp.clear();
            for (i=0;i<4;i++) bd[i]=data.a[Mindex.at(P-1)].bd[i]-bd[i];
            if (bond==1) {bd[1]+=1;bd[0]-=1;}
            else if (bond==2) {bd[2]+=1;bd[0]-=1;}
            k=0;j=0;i=id2;
            if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[Mindex[P-1]].chg==data.a[i].chg) {Mindex.at(P-1)=i;Rindex.at(n)=R+bond;}
            else {
				Rindex.at(n)=R;
                Mindex.at(n)=M;
                reset();
                return 0;
            }
        }
    reset();
    return 1;
}
int MOLECULE::chk_chg(MOLECULE &mol) {
	rechg();
	mol.rechg();
	int num[2];
	char dot='.';
	ionic="";
	num[0]=abs(mol.chg);
	num[1]=abs(chg);
	if (num[0]==0 && num[1]==0) ionic=molesmi;
	else if (num[0]==0) ionic=mol.molesmi;
	else if (num[1]==0) ionic=molesmi;
	else {
		int i;
		for (i=0;i<num[0];i++) {
			ionic+="(";
			ionic+=molesmi;
			ionic+=")";
		}
		ionic+=dot;
		for (i=0;i<num[1];i++) {
			ionic+="(";
            ionic+=mol.molesmi;
			ionic+=")";
        }
	}
	return 1;
}
int MOLECULE::rechg() {
	int i=0,j=0;
	chg=0;
	for (j=0;j<Cindex.size();j++) {
		chg+=data.a[Mindex.at(j)].chg;
	}
	return 1;
}

void MOLECULE::neutralize() {
	int i,j,k,m,n;
	chg=0;
	for (i=0;i<Cindex.size();i++) {
		chg+=data.a[Mindex.at(i)].chg;
		if (data.a[Mindex.at(i)].chg) j=i;
	}
	if (chg) {
		if (chg<0) {
			m=0;
			while (1) {
				m++;
				n=0;
				k=rand()%Cindex.size();
				for (i=0;i<data.a[Mindex.at(k)].nbond;i++) {
					if (Bindex[k][i] && Bindex[k][i]<3) n=Bindex[k][i];
				}
				if (k!=j && n) break;
				if (m>5) return;
			}
			if (n==2) {
				if (rand()%2) {
					Cindex.push_back(Mindex.size());
					Pindex.push_back(Cindex.at(k));
					Mindex.push_back(16);
					Cyindex.push_back(0);
					Rindex.push_back(2);
				}
				else {
                    Cindex.push_back(Mindex.size());
                    Pindex.push_back(Cindex.at(k));
                    Mindex.push_back(18);
                    Cyindex.push_back(0);
                    Rindex.push_back(2);
				}
			}
			else if (n==1){
				Cindex.push_back(Mindex.size());
				Pindex.push_back(Cindex.at(k));
				Mindex.push_back(rand()%3+15);
				Cyindex.push_back(0);
				Rindex.push_back(1);
			}
		}
		else if (chg>0) {
			m=0;
			while (1) {
                n=0;
				m++;
                k=rand()%Cindex.size();
                for (i=0;i<data.a[Mindex.at(k)].nbond;i++) {
                    if (Bindex[k][i] && Bindex[k][i]==1) n=Bindex[k][i];
                }
                if (k!=j && n) break;
				if (m>5) return;
            }
			if (n) {
				int tmp[3]={30,29,35};
				Cindex.push_back(Mindex.size());
				Pindex.push_back(Cindex.at(k));
				Mindex.push_back(tmp[rand()%3]);
				Cyindex.push_back(0);
				Rindex.push_back(1);
			}
		}
	}
		reset();
		mds2smi();
		return;
}

void MOLECULE::mds23d(ostream &outs) {
	int i,j,k,n;
	double sum;
	Matrix vec[Mindex.size()][6];
	int **cont;
	cont=new int *[Mindex.size()];
	for (i=0;i<Mindex.size();i++) cont[i]=new int [6];
	vector<double> x[3];
	vector<int> mm;
	vector<OPT> opt;
	opt.resize(1);
	for (i=0;i<Mindex.size();i++) mm.push_back(Mindex.at(i));
	for (i=0;i<mm.size();i++) {
		if (mm.at(i)==10) mm.at(i)=5;
		else if (mm.at(i)==28) mm.at(i)=29;
	}
	for (i=0;i<mm.size();i++) {
		k=data.a[mm.at(i)].norder;
		n=0;
		for (j=0;j<data.a[mm.at(i)].norder;j++) {
			if (Bindex[i][j]>1) {
				k+=(Bindex[i][j]-1);
				n=1;
			}
		}
		if (k>data.a[Mindex.at(i)].norder) {
			for (j=1;j<data.num;j++) {
				if (n==1 && data.a.at(j).atm==data.a.at(mm.at(i)).atm && data.a.at(j).norder==k && data.a.at(j).chg==data.a.at(mm.at(i)).chg) {
					mm.at(i)=j;
					break;
				}
			}
		}
	}
	for (i=0;i<Mindex.size();i++) {
		for (j=0;j<6;j++) {
			vec[i][j].resize(1,1);
			cont[i][j]=0;
		}
	}
	
	for (i=0;i<mm.size();i++) {
		for (j=0;j<data.a[mm.at(i)].norder;j++) {
			vec[i][j].resize(3,1);
			cont[i][j]=1;
		}
	}
	/* initailize */
	if (data.a[mm.at(0)].type==1) {
		vec[0][0].Mtrx[0][0]=0.000000;
		vec[0][0].Mtrx[1][0]=0.000000;
		vec[0][0].Mtrx[2][0]=1.000000;
		vec[0][1].Mtrx[0][0]=-0.78335;
		vec[0][1].Mtrx[1][0]=-0.4523;
		vec[0][1].Mtrx[2][0]=-0.4264;
		vec[0][2].Mtrx[0][0]=0.000000;
		vec[0][2].Mtrx[1][0]=0.9045;
		vec[0][2].Mtrx[2][0]=-0.4264;
		vec[0][3].Mtrx[0][0]=0.78335;
		vec[0][3].Mtrx[1][0]=-0.4523;
		vec[0][3].Mtrx[2][0]=-0.4264;
		vec[0][0].rotate(vec[0][1],vec[0][0],15.0*3.1415926/180.0);
		vec[0][3].rotate(vec[0][1],vec[0][3],15.0*3.1415926/180.0);
		vec[0][2].rotate(vec[0][1],vec[0][2],15.0*3.1415926/180.0);
	}
	else if (data.a[mm.at(0)].type==2) {
		Matrix ref(3,1);
		ref[0][0]=ref[1][0]=ref[2][0]=1.0;
		vec[0][0].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[2][0]=1.0/sqrt(3.0);
		vec[0][1].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][1].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][1].Mtrx[2][0]=1.0/sqrt(3.0);
		vec[0][2].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][2].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][2].Mtrx[2][0]=1.0/sqrt(3.0);
		vec[0][1].rotate(ref,vec[0][1],120.0*3.1415926/180.0);
		vec[0][2].rotate(ref,vec[0][2],-120.0*3.1415926/180.0);
	}
	else if (data.a[mm.at(0)].type==3) {
        Matrix ref(3,1);
        ref[0][0]=ref[1][0]=ref[2][0]=1.0;
        vec[0][0].Mtrx[0][0]=1.0/sqrt(3.0);
        vec[0][0].Mtrx[1][0]=1.0/sqrt(3.0);
        vec[0][0].Mtrx[2][0]=1.0/sqrt(3.0);
        vec[0][1].Mtrx[0][0]=1.0/sqrt(3.0);
        vec[0][1].Mtrx[1][0]=1.0/sqrt(3.0);
        vec[0][1].Mtrx[2][0]=1.0/sqrt(3.0);
        vec[0][2].Mtrx[0][0]=1.0/sqrt(3.0);
        vec[0][2].Mtrx[1][0]=1.0/sqrt(3.0);
        vec[0][2].Mtrx[2][0]=1.0/sqrt(3.0);
        vec[0][1].rotate(ref,vec[0][1],120.0*3.1415926/180.0);
        vec[0][2].rotate(ref,vec[0][2],-120.0*3.1415926/180.0);

	}
	else if (data.a[mm.at(0)].type==4) {
		vec[0][0].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[2][0]=1.0/sqrt(3.0);
		vec[0][1].Mtrx[0][0]=-1.0/2.0;
		vec[0][1].Mtrx[1][0]=1.0/2.0;
		vec[0][1].Mtrx[2][0]=-1.0/sqrt(2.0);
	}
	else if (data.a[mm.at(0)].type==5) {
		vec[0][0].Mtrx[0][0]=0.3;
		vec[0][0].Mtrx[1][0]=0.1;
		vec[0][0].Mtrx[2][0]=1.0000;
		vec[0][1].Mtrx[0][0]=0.3;
		vec[0][1].Mtrx[1][0]=0.1;
		vec[0][1].Mtrx[2][0]=-1.0000;
        vec[0][2].Mtrx[0][0]=-0.50000*sqrt(3.0000);
        vec[0][2].Mtrx[1][0]=-sqrt(3.0000)/6.0000*sqrt(3.00000);
        vec[0][2].Mtrx[2][0]=0.100000;
        vec[0][3].Mtrx[0][0]=0.50000*sqrt(3.0000);
        vec[0][3].Mtrx[1][0]=-sqrt(3.0000)/6.0000*sqrt(3.00000);
        vec[0][3].Mtrx[2][0]=0.300000;
        vec[0][4].Mtrx[0][0]=0.1000000;
        vec[0][4].Mtrx[1][0]=1.0000;
        vec[0][4].Mtrx[2][0]=0.30000;
		for (i=0;i<5;i++) {
			sum=0.0;
            for (j=0;j<3;j++) sum+=vec[0][i].Mtrx[j][0]*vec[0][i].Mtrx[j][0];
            sum=sqrt(sum);
            for (j=0;j<3;j++) vec[0][i].Mtrx[j][0]=vec[0][i].Mtrx[j][0]/sum;
		}
	}
	else if (data.a[mm.at(0)].type==6) {
        vec[0][0].Mtrx[0][0]=0.2;
        vec[0][0].Mtrx[1][0]=0.1;
        vec[0][0].Mtrx[2][0]=1.0000;
        vec[0][1].Mtrx[0][0]=0.3;
        vec[0][1].Mtrx[1][0]=0.5;
        vec[0][1].Mtrx[2][0]=-1.0000;
		vec[0][2].Mtrx[0][0]=1.00000;
		vec[0][2].Mtrx[1][0]=0.1;
		vec[0][2].Mtrx[2][0]=0.3;
		vec[0][3].Mtrx[0][0]=-1.0000;
		vec[0][3].Mtrx[1][0]=0.5;
		vec[0][3].Mtrx[2][0]=0.3;
		vec[0][4].Mtrx[0][0]=0.3;
		vec[0][4].Mtrx[1][0]=1.000000;
		vec[0][4].Mtrx[2][0]=0.3;
		vec[0][5].Mtrx[0][0]=0.2;
		vec[0][5].Mtrx[1][0]=-1.0000;
		vec[0][5].Mtrx[2][0]=0.3;
		for (i=0;i<6;i++) {
			sum=0.0;
			for (j=0;j<3;j++) sum+=vec[0][i].Mtrx[j][0]*vec[0][i].Mtrx[j][0];
			sum=sqrt(sum);
			for (j=0;j<3;j++) vec[0][i].Mtrx[j][0]=vec[0][i].Mtrx[j][0]/sum;
		}
	}
	else if (data.a[mm.at(0)].type==7) {
		vec[0][0].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[2][0]=1.0/sqrt(3.0);
	}
	x[0].push_back(0.0000000000);
	x[1].push_back(0.0000000000);
	x[2].push_back(0.0000000000);
	int P=0;
	double *pt;
	pt=new double [3];
    double *xx;
	xx=new double [3];
	int memo;
	for (i=1;i<mm.size();i++) {
		P=Pindex.at(i)-1;
		for (j=0;j<data.a[mm.at(P)].norder;j++) {
			if (cont[P][j]) {
				cont[P][j]=0;
				memo=j;
				break;
			}
		}
		cont[i][0]=0;
		for (k=0;k<3;k++) x[k].push_back(x[k].at(P)+(data.a[mm.at(P)].rb+data.a[mm.at(i)].rb)*vec[P][j].Mtrx[k][0]);
		for (k=0;k<3;k++) vec[i][0].Mtrx[k][0]=-vec[P][j].Mtrx[k][0];
		if (data.a[mm.at(i)].type==1) { /* tetrahedral */
			for (k=0;k<3;k++) pt[k]=-0.5000*vec[i][0].Mtrx[k][0];
			int tmpp[3];
			for (k=0;k<3;k++) tmpp[3]=(fabs(vec[i][0].Mtrx[k][0])>1E-5);
			for (k=0;k<3;k++) {
				if (tmpp[k]) xx[k]=pt[k]+(double)rand()/RAND_MAX*3.0;
				else xx[k]=pt[k];
			}
			sum=0.0;
			if (tmpp[0]) {
				sum=vec[i][0].Mtrx[1][0]*(pt[1]-xx[1])+vec[i][0].Mtrx[2][0]*(pt[2]-xx[2]);
				xx[0]=sum/vec[i][0].Mtrx[0][0]+pt[0];
			}
			else if (tmpp[1]) {
				sum=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])+vec[i][0].Mtrx[2][0]*(pt[2]-xx[2]);
                xx[1]=sum/vec[i][0].Mtrx[1][0]+pt[1];
			}
			else {
                sum=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])+vec[i][0].Mtrx[1][0]*(pt[1]-xx[1]);
                xx[1]=sum/vec[i][0].Mtrx[2][0]+pt[2];
			}	
			sum=0;
			sum=(xx[0]-pt[0])*(xx[0]-pt[0])+(xx[1]-pt[1])*(xx[1]-pt[1])+(xx[2]-pt[2])*(xx[2]-pt[2]);
			sum=sqrt(sum);
			sum/=(sqrt(3.0)/2.0);
			xx[0]=(xx[0]-pt[0])/sum;
			xx[1]=(xx[1]-pt[1])/sum;
			xx[2]=(xx[2]-pt[2])/sum;
            for (k=1;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].Mtrx[0][0]=xx[0];
				vec[i][k].Mtrx[1][0]=xx[1];
				vec[i][k].Mtrx[2][0]=xx[2];
			}
			for (k=2;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].rotate(vec[i][0],vec[i][k],2.0000*3.14159/3.0000*(double)pow(-1,k));
				sum=0;
				sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
				sum=sqrt(sum);
				vec[i][k].Mtrx[0][0]/=sum;
				vec[i][k].Mtrx[1][0]/=sum;
				vec[i][k].Mtrx[2][0]/=sum;
			}
		}
		else if (0) {
            for (k=0;k<3;k++) pt[k]=x[k].at(i)-0.5000*vec[i][0].Mtrx[k][0];
            for (k=0;k<2;k++) {
                xx[k]=pt[k]+(double)rand()/RAND_MAX*3.0;
            }
            sum=0;
            sum=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])+vec[i][0].Mtrx[1][0]*(pt[1]-xx[1]);
            if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[2]=sum/vec[i][0].Mtrx[2][0]+pt[2];
            else {
				xx[2]=pt[2];
                if (fabs(vec[i][0].Mtrx[1][0])>1E-5) xx[1]=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])/vec[i][0].Mtrx[1][0]+pt[1];
                else xx[1]=pt[1];
                if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[0]=vec[i][0].Mtrx[1][0]*(pt[1]-xx[1])/vec[i][0].Mtrx[0][0]+pt[0];
                else xx[0]=pt[0];
			}
            sum=0;
            sum=(xx[0]-pt[0])*(xx[0]-pt[0])+(xx[1]-pt[1])*(xx[1]-pt[1])+(xx[2]-pt[2])*(xx[2]-pt[2]);
            sum=sqrt(sum);
            sum=sum/(sqrt(3.0)/2.0);
            xx[0]=(xx[0]-pt[0])/sum;
            xx[1]=(xx[1]-pt[1])/sum;
            xx[2]=(xx[2]-pt[2])/sum;
            for (k=1;k<data.a[mm.at(i)].norder;k++) {
                vec[i][k].Mtrx[0][0]=xx[0]-0.50000*vec[i][0].Mtrx[0][0];
                vec[i][k].Mtrx[1][0]=xx[1]-0.50000*vec[i][0].Mtrx[1][0];
                vec[i][k].Mtrx[2][0]=xx[2]-0.50000*vec[i][0].Mtrx[2][0];
            }
            for (k=2;k<data.a[mm.at(i)].norder;k++) {
                vec[i][k].rotate(vec[i][0],vec[i][k],2.0000*3.14159/3.0000*(double)pow(-1,k));
                sum=0;
                sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
                sum=sqrt(sum);
                vec[i][k].Mtrx[0][0]/=sum;
                vec[i][k].Mtrx[1][0]/=sum;
                vec[i][k].Mtrx[2][0]/=sum;
            }
        }
		else if (data.a[mm.at(i)].type==3||data.a[mm.at(i)].type==2) {
			Matrix a(3,1),aa(3,1);
			if (data.a[mm.at(P)].norder>2) {
        		for (j=0;j<data.a[mm.at(P)].norder;j++) {
            		if (j!=memo) break;
            	}
				n=0;
				a[0][0]=vec[i][0].Mtrx[1][0]*vec[P][j].Mtrx[2][0]-vec[i][0].Mtrx[2][0]*vec[P][j].Mtrx[1][0];
				a[1][0]=vec[i][0].Mtrx[2][0]*vec[P][j].Mtrx[1][0]-vec[i][0].Mtrx[1][0]*vec[P][j].Mtrx[2][0];
				a[2][0]=vec[i][0].Mtrx[0][0]*vec[P][j].Mtrx[1][0]-vec[i][0].Mtrx[1][0]*vec[P][j].Mtrx[0][0];
				if (fabs(a[0][0])<1E-3 && fabs(a[1][0])<1E-3 && fabs(a[2][0])<1E-3) n=1;
				if (n) {
					while (1) {
                    	for (k=0;k<3;k++) aa[k][0]=(double)rand()/RAND_MAX;
                    	a[0][0]=aa[1][0]*vec[P][j].Mtrx[2][0]-aa[2][0]*vec[P][j].Mtrx[1][0];
                    	a[1][0]=aa[2][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[2][0];
                    	a[2][0]=aa[0][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[0][0];
                    	if (fabs(a[0][0])>1E-3 || fabs(a[1][0])>1E-3 || fabs(a[2][0])>1E-3) break;
                	}
				}
				for (k=1;k<data.a[mm.at(i)].norder;k++) {
                	vec[i][k].Mtrx[0][0]=vec[i][0].Mtrx[0][0];
                	vec[i][k].Mtrx[1][0]=vec[i][0].Mtrx[1][0];
                	vec[i][k].Mtrx[2][0]=vec[i][0].Mtrx[2][0];
            	}
            	for (k=1;k<data.a[mm.at(i)].norder;k++) {
                	vec[i][k].rotate(a,vec[i][k],2.0000*3.14159/3.0000*(double)pow(-1,k));
                	sum=0;
                	sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
                	sum=sqrt(sum);
                	vec[i][k].Mtrx[0][0]/=sum;
                	vec[i][k].Mtrx[1][0]/=sum;
                	vec[i][k].Mtrx[2][0]/=sum;
            	}
			}
			else {
				while (1) {
					for (k=0;k<3;k++) aa[k][0]=(double)rand()/RAND_MAX;
					a[0][0]=aa[1][0]*vec[P][j].Mtrx[2][0]-aa[2][0]*vec[P][j].Mtrx[1][0];
                	a[1][0]=aa[2][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[2][0];
                	a[2][0]=aa[0][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[0][0];
					if (fabs(a[0][0])>1E-3 || fabs(a[1][0])>1E-3 || fabs(a[2][0])>1E-3) break;
				}
				for (k=1;k<data.a[mm.at(i)].norder;k++) {
                    vec[i][k].Mtrx[0][0]=vec[i][0].Mtrx[0][0];
                    vec[i][k].Mtrx[1][0]=vec[i][0].Mtrx[1][0];
                    vec[i][k].Mtrx[2][0]=vec[i][0].Mtrx[2][0];
                }
                for (k=1;k<data.a[mm.at(i)].norder;k++) {
                    vec[i][k].rotate(a,vec[i][k],2.0000*3.14159/3.0000*(double)pow(-1,k));
                    sum=0;
                    sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
                    sum=sqrt(sum);
                    vec[i][k].Mtrx[0][0]/=sum;
                    vec[i][k].Mtrx[1][0]/=sum;
                    vec[i][k].Mtrx[2][0]/=sum;
				}
			}
        }
		else if (data.a[mm.at(i)].type==4) {
			Matrix a(3,1);
			Matrix aa(3,1);
			while (1) {
                    for (k=0;k<3;k++) aa[k][0]=(double)rand()/RAND_MAX;
                    a[0][0]=aa[1][0]*vec[P][j].Mtrx[2][0]-aa[2][0]*vec[P][j].Mtrx[1][0];
                    a[1][0]=aa[2][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[2][0];
                    a[2][0]=aa[0][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[0][0];
                    if (fabs(a[0][0])>1E-3 || fabs(a[1][0])>1E-3 || fabs(a[2][0])>1E-3) break;
          	}
			for (k=1;k<data.a[mm.at(i)].norder;k++) {
                    vec[i][k].Mtrx[0][0]=vec[i][0].Mtrx[0][0];
                    vec[i][k].Mtrx[1][0]=vec[i][0].Mtrx[1][0];
                    vec[i][k].Mtrx[2][0]=vec[i][0].Mtrx[2][0];
            }
			for (k=1;k<data.a[mm.at(i)].norder;k++) {
                    vec[i][k].rotate(a,vec[i][k],2.0000*3.14159/3.0000*(double)pow(-1,k));
                    sum=0;
                    sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
                    sum=sqrt(sum);
                    vec[i][k].Mtrx[0][0]/=sum;
                    vec[i][k].Mtrx[1][0]/=sum;
                    vec[i][k].Mtrx[2][0]/=sum;
            }

		}
        else if (data.a[mm.at(i)].type==5) {
            for (k=0;k<3;k++) pt[k]=x[k].at(i);
            for (k=0;k<2;k++) {
                xx[k]=pt[k]+(double)rand()/RAND_MAX*3.0;
            }
            sum=0;
            sum=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])+vec[i][0].Mtrx[1][0]*(pt[1]-xx[1]);
            if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[2]=sum/vec[i][0].Mtrx[2][0]+pt[2];
            else {
                xx[2]=pt[2];
                if (fabs(vec[i][0].Mtrx[1][0])>1E-5) xx[1]=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])/vec[i][0].Mtrx[1][0]+pt[1];
                else xx[1]=pt[1];
                if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[0]=vec[i][0].Mtrx[1][0]*(pt[1]-xx[1])/vec[i][0].Mtrx[0][0]+pt[0];
                else xx[0]=pt[0];
            }
            sum=0;
            sum=(xx[0]-pt[0])*(xx[0]-pt[0])+(xx[1]-pt[1])*(xx[1]-pt[1])+(xx[2]-pt[2])*(xx[2]-pt[2]);
            sum=sqrt(sum);
            xx[0]=(xx[0]-pt[0])/sum;
            xx[1]=(xx[1]-pt[1])/sum;
            xx[2]=(xx[2]-pt[2])/sum;
            for (k=1;k<data.a[mm.at(i)].norder;k++) {
                vec[i][k].Mtrx[0][0]=xx[0];
                vec[i][k].Mtrx[1][0]=xx[1];
                vec[i][k].Mtrx[2][0]=xx[2];
            }
            for (k=2;k<data.a[mm.at(i)].norder;k++) {
                vec[i][k].rotate(vec[i][0],vec[i][k],2.0000*3.1415926/4.00000*((double)k-1.00000));
                sum=0;
                sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
                sum=sqrt(sum);
                vec[i][k].Mtrx[0][0]/=sum;
                vec[i][k].Mtrx[1][0]/=sum;
                vec[i][k].Mtrx[2][0]/=sum;
            }
        }
		else if (data.a[mm.at(i)].type==6) {
            for (k=0;k<3;k++) pt[k]=x[k].at(i);
            for (k=0;k<2;k++) {
                xx[k]=pt[k]+(double)rand()/RAND_MAX*3.0;
            }
            sum=0;
            sum=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])+vec[i][0].Mtrx[1][0]*(pt[1]-xx[1]);
            if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[2]=sum/vec[i][0].Mtrx[2][0]+pt[2];
            else {
                xx[2]=pt[2];
                if (fabs(vec[i][0].Mtrx[1][0])>1E-5) xx[1]=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])/vec[i][0].Mtrx[1][0]+pt[1];
                else xx[1]=pt[1];
                if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[0]=vec[i][0].Mtrx[1][0]*(pt[1]-xx[1])/vec[i][0].Mtrx[0][0]+pt[0];
                else xx[0]=pt[0];
            }
            sum=0;
            sum=(xx[0]-pt[0])*(xx[0]-pt[0])+(xx[1]-pt[1])*(xx[1]-pt[1])+(xx[2]-pt[2])*(xx[2]-pt[2]);
            sum=sqrt(sum);
            xx[0]=(xx[0]-pt[0])/sum;
            xx[1]=(xx[1]-pt[1])/sum;
            xx[2]=(xx[2]-pt[2])/sum;
            for (k=1;k<data.a[mm.at(i)].norder;k++) {
                vec[i][k].Mtrx[0][0]=xx[0];
                vec[i][k].Mtrx[1][0]=xx[1];
                vec[i][k].Mtrx[2][0]=xx[2];
            }
            for (k=2;k<data.a[mm.at(i)].norder;k++) {
                vec[i][k].rotate(vec[i][0],vec[i][k],108.5*3.1415926/180.0*((double)k-1.000));
                sum=0;
                sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
                sum=sqrt(sum);
                vec[i][k].Mtrx[0][0]/=sum;
                vec[i][k].Mtrx[1][0]/=sum;
                vec[i][k].Mtrx[2][0]/=sum;
            }
		}
		else if (data.a[mm.at(i)].type==7) {}
			
	}
	int q=0;
	for (i=0;i<mm.size();i++) {
		if (Cyindex.at(i)) {
			for (j=0;j<data.a[mm.at(i)].norder;j++) {
				if (cont[i][j]) {
					for (k=i+1;k<mm.size();k++) {
						if (Cyindex.at(i)==Cyindex.at(k)) {
							for (n=0;n<data.a[mm.at(k)].norder;n++) {
								if (cont[k][n]) {
									sum=0.0;
									sum=fabs(vec[i][j].Mtrx[0][0]+vec[k][n].Mtrx[0][0])+fabs(vec[i][j].Mtrx[1][0]+vec[k][n].Mtrx[1][0])+fabs(vec[i][j].Mtrx[2][0]+vec[k][n].Mtrx[2][0]);
									if (sum<1E-3) {
										cont[i][j]=cont[k][n]=0;
										q=1;
										break;
									}
								}
							}
						}
					}
					if (q) break;
				}
			}
			if (q) break;
		}
		if (q) break;
	}
	if (!q) {
		for (i=0;i<mm.size();i++) {
			if (Cyindex.at(i)) {
				for (j=0;j<data.a[mm.at(i)].norder;j++) {
					if (cont[i][j]) {
						cont[i][j]=0;
						break;
					}
				}
			}
		}
	}
	for (i=0;i<mm.size();i++) {
		for (j=0;j<data.a[mm.at(i)].norder;j++) {
			if (cont[i][j]) {
				for (k=0;k<3;k++) {
					x[k].push_back(x[k].at(i)+vec[i][j].Mtrx[k][0]*(data.a[mm.at(i)].rb+0.32));
				}
				cont[i][j]=0;
			}
		}
	}
	opt[0].count=if_circle;
	opt[0].num=x[0].size();
	opt[0].f.resize(opt[0].num*3,1);
	opt[0].df.resize(1,opt[0].num*3);
	opt[0].x0.resize(opt[0].num*3,1);
	opt[0].id=new int [opt[0].num];
	opt[0].table=new int *[opt[0].num];
	for (i=0;i<opt[0].num;i++) opt[0].table[i]=new int [opt[0].num];
	for (i=0;i<opt[0].num;i++) {
		if (i<mm.size()) opt[0].id[i]=mm.at(i);
		else opt[0].id[i]=0;
	}
	for (i=0;i<opt[0].num;i++) for (j=0;j<opt[0].num;j++) opt[0].table[i][j]=0;
	for (i=1;i<mm.size();i++) {
		opt[0].table[i][Pindex.at(i)-1]=opt[0].table[Pindex.at(i)-1][i]=1;
	}
    if (if_circle) {
        int c1[2],c2[2];
        for (i=0;i<Cyindex.size();i++) {
            if (Cyindex.at(i)) {
                if (Cyindex.at(i)>10) {
                    c1[0]=Cyindex.at(i)/10;
                    c1[1]=Cyindex.at(i)%10;
                }
                else c1[0]=c1[1]=Cyindex.at(i);
                for (j=i+1;j<Cyindex.size();j++) {
                    if (Cyindex.at(j)>10) {
                        c2[0]=Cyindex.at(j)/10;
                        c2[1]=Cyindex.at(j)%10;
                    }
                    else c2[0]=c2[1]=Cyindex.at(j);
                    if (c1[0]==c2[0] || c1[0]==c2[1] || c1[1]==c2[0] || c1[1]==c2[1]) {
                        opt[0].table[i][j]=opt[0].table[j][i]=1;
                    }
                }
            }
        }
    }

	int t=0;
	for (i=0;i<mm.size();i++) {
		k=0;
		for (j=0;j<opt[0].num;j++) k+=opt[0].table[i][j];
		if (k<data.a[mm.at(i)].norder) {
			for (n=0;n<(data.a[mm.at(i)].norder-k);n++) {
				opt[0].table[i][mm.size()+t]=opt[0].table[mm.size()+t][i]=1;
				t++;
			}
		}
	}
	for (i=0;i<opt[0].num;i++) {
		for (j=0;j<3;j++) {
			opt[0].x0[i+j*opt[0].num][0]=x[j].at(i);
		}
	}	
		opt[0].optimize();
		opt[0].output(outs);
		outs<<endl;
		delete [] opt[0].id;
		for (i=0;i<opt[0].num;i++) {
			delete [] opt[0].table[i];
		}
		delete [] opt[0].table;
		for (i=0;i<Mindex.size();i++) delete [] cont[i];
		delete [] cont;
		vector<int>().swap(mm);
		vector<double>().swap(x[0]);
		vector<double>().swap(x[1]);
		vector<double>().swap(x[2]);
		vector<OPT>().swap(opt);
		delete [] pt;
		delete [] xx;
	return;
}
