#include "ATOM.h"
using namespace std;

void DEATOM::find_r() {
	if (name=="C") {
		r_bnd=0.77;
		nbnd=4;
	}
	else if (name=="N") {
		r_bnd=0.74;
		nbnd=3;
	}
	else if (name=="O") {
		r_bnd=0.74;
		nbnd=2;
	}
	else if (name=="F") {
		r_bnd=0.72;
		nbnd=1;
	}
	else if (name=="Cl") {
		r_bnd=0.99;
		nbnd=1;;
	}
	else if (name=="Br") {
		r_bnd=1.20; 
		nbnd=1;
	}
	else if (name=="I") {
		r_bnd=1.39;
		nbnd=1;
	}
	else if (name=="H") {
		r_bnd=0.37;
		nbnd=1;
	}
	else if (name=="S") {
		r_bnd=1.05;
		nbnd=2;
	}
	else if (name=="P") {
		r_bnd=1.07;
		nbnd=4;
	}
	
	return;
}

void POOL::set_up() {
	num=36;
	int i,j,k;
	a.resize(num);
	for (i=0;i<num;i++) a.at(i).order.resize(6);
	for (i=0;i<num;i++) {
			a.at(i).nh=0;
			for (j=0;j<6;j++) a.at(i).order.at(j)=0;
			for (j=0;j<3;j++) a.at(i).bd[j]=0;
			a.at(i).index=2; // for neutral atom
			a.at(i).chg=0;	// for neutral atom
	}
	double f=3.1415926/180.0;
	a[0].id=0;
	a[0].name="H(-)";
	a[0].nbond=4;
	a[0].order[0]=1;
	a[0].order[1]=0;
	a[0].order[2]=0;
	a[0].order[3]=0;
	a[0].bd[0]=1;
	a[0].nbond=1;
	a[0].norder=1;
	a[0].rb=0.37;
	a[0].ang0=0.0;
	a[0].atm="H";

	a[1].id=1;
	a[1].name="C(-)(-)(-)(-)";
	a[1].order[0]=a[1].order[1]=a[1].order[2]=a[1].order[3]=1;
	a[1].nbond=13;
	a[1].norder=4;
	a[1].rb=0.77;
	a[1].type=1;
	a[1].bd[0]=4;
	a[1].ang0=109.5*f;
	a[1].atm="C";

	a[2].id=2;
	a[2].name="C(=)(-)(-)";
	a[2].order[0]=2;
	a[2].order[1]=a[2].order[2]=1;
	a[2].order[3]=0;
	a[2].nbond=10;
	a[2].norder=3;
	a[2].rb=0.77;
	a[2].type=3;
	a[2].bd[0]=2;
	a[2].bd[1]=1;
	a[2].ang0=120.0*f;
	a[2].atm="C";

	a[3].id=3;
	a[3].name="C(#)(-)";
	a[3].order[0]=3;
	a[3].order[1]=1;
	a[3].order[2]=a[3].order[3]=0;
	a[3].nbond=7;
	a[3].norder=2;
	a[3].rb=0.77;
	a[3].type=4;
	a[3].bd[0]=1;
	a[3].bd[2]=1;
	a[3].ang0=180*f;
	a[3].atm="C";

	a[4].id=4;
	a[4].name="C(=)(=)";
	a[4].order[0]=a[4].order[1]=2;
	a[4].order[2]=a[4].order[3]=0;
	a[4].nbond=7;
	a[4].norder=2;
	a[4].rb=0.77;
	a[4].type=4;
	a[4].bd[1]=2;
	a[4].ang0=180*f;
	a[4].atm="C";

	a[5].id=5;
	a[5].name="O(-)(-)";
	a[5].order[0]=a[5].order[1]=1;
	a[5].order[2]=a[5].order[3]=0;
	a[5].nbond=7;
	a[5].norder=2;
	a[5].rb=0.74;
	a[5].type=4;
	a[5].bd[0]=2;
	a[5].ang0=104.5*f;
	a[5].atm="O";

	a[6].id=6;
	a[6].name="O(=)";
	a[6].order[0]=2;
	a[6].order[1]=a[6].order[2]=a[6].order[3]=0;
	a[6].nbond=4;
	a[6].norder=1;
	a[6].rb=0.74;
	a[6].type=7;
	a[6].bd[1]=1;
	a[6].ang0=0;
	a[6].atm="O";

	a[7].id=7;
	a[7].name="N(-)(-)(-)";
	a[7].order[0]=a[7].order[1]=a[7].order[2]=1;
	a[7].order[3]=0;
	a[7].nbond=10;
	a[7].norder=3;
	a[7].rb=0.74;
	a[7].type=2;
	a[7].bd[0]=3;
	a[7].ang0=107.5*f;
	a[7].atm="N";

	a[8].id=8;
	a[8].name="N(=)(-)";
	a[8].order[0]=2;
	a[8].order[1]=1;
	a[8].order[2]=a[8].order[3]=0;
	a[8].nbond=7;
	a[8].norder=2;
	a[8].rb=0.74;
	a[8].type=4;
	a[8].bd[0]=1;
	a[8].bd[1]=1;
	a[8].ang0=115.0*f;
	a[8].atm="N";

	a[9].id=9;
	a[9].name="N(#)";
	a[9].order[0]=3;
	a[9].order[1]=a[9].order[2]=a[9].order[3]=0;
	a[9].nbond=4;
	a[9].norder=1;
	a[9].rb=0.74;
	a[9].type=7;
	a[9].bd[2]=1;
	a[9].ang0=0;
	a[9].atm="N";

	a[10].id=10;
	a[10].name="O(-)";
	a[10].order[0] = 1;
	a[10].order[1] = a[10].order[2] =a[10].order[3] =0;
	a[10].nbond=4;
	a[10].norder=1;
	a[10].nh=1;
	a[10].rb=0.74;
	a[10].type=4;
	a[10].bd[0]=1;
	a[10].ang0=0;
	a[10].atm="O";

	a[11].id = 11;
	a[11].name = "F(-)";
	a[11].order[0] = 1;
	a[11].order[1] = a[11].order[2] =a[11].order[3] =0;
	a[11].nbond=4;
	a[11].norder=1;
	a[11].rb=0.72;
	a[11].type=7;
	a[11].bd[0]=1;
	a[11].ang0=0;
	a[11].atm="F";

	a[12].id = 12;
	a[12].name="Cl(-)";
	a[12].order[0] = 1;
	a[12].order[1]=a[12].order[2]=a[12].order[3] = 0;
	a[12].nbond=5;
	a[12].norder=1;
	a[12].atm="Cl";	

	a[13].id = 13;
    a[13].name="Br(-)";
    a[13].order[0] = 1;
    a[13].order[1] = a[13].order[2] = a[13].order[3] =0;
    a[13].nbond = 5;
    a[13].norder=1;
	a[13].rb=1.200000;
	a[13].type=7;
	a[13].bd[0]=1;
	a[13].index=3;
	a[13].ang0=0;
	a[13].atm="Br";

	a[14].id = 14;
	a[14].name ="I(-)";
	a[14].order[0]=1;
	a[14].order[1] = a[14].order[2] = a[14].order[3] = 0;
	a[14].nbond = 4;
	a[14].norder=1;
	a[14].rb=1.3900000;
	a[14].type=7;
	a[14].bd[0]=1;
	a[14].ang0=0;
	a[14].atm="I";
	
	a[15].id=15;
	a[15].name="[NH0+](-)(-)(-)(-)";
	for (i=0;i<4;i++) a[15].order[i]=1;
	a[15].nbond=18;
	a[15].norder=4;
	a[15].rb=0.74;
	a[15].type=1;
	a[15].bd[0]=4;
	a[15].chg=1;
	a[15].index=7;
	a[15].ang0=109.5*f;
	a[15].atm="N";

	a[16].id=16;
	a[16].name="[NH0+](=)(-)(-)";
	a[16].order[0]=2;
	a[16].order[1]=a[16].order[2]=1;
	a[16].order[3]=0;
	a[16].nbond=15;
	a[16].norder=3;
	a[16].rb=0.74;
	a[16].type=3;
	a[16].bd[0]=2;
	a[16].bd[1]=1;
	a[16].chg=1;
	a[16].index=7;
	a[16].ang0=120.0*f;
	a[16].atm="N";

	a[17].id=17;
	a[17].name="[PH0+](-)(-)(-)(-)";
	for (i=0;i<4;i++) a[17].order[i]=1;
	a[17].norder=4;
	a[17].nbond=18;
	a[17].rb=1.07;
	a[17].type=1;
	a[17].bd[0]=4;
	a[17].chg=1;
	a[17].index=7;
	a[17].ang0=109.5*f;
	a[17].atm="P";

	a[18].id=18;
	a[18].name="[PH0+](=)(-)(-)";
	a[18].order[0]=2;
	a[18].order[1]=a[18].order[2]=1;
	a[18].order[3]=0;
	a[18].norder=3;
	a[18].nbond=15;
	a[18].rb=1.07;
	a[18].bd[0]=2;
	a[18].bd[1]=1;
	a[18].chg=1;
	a[18].index=7;
	a[18].ang0=120.0*f;
	a[18].atm="P";

	a[19].id=19;
	a[19].name="S(-)(-)";
	a[19].order[0]=a[19].order[1]=1;
	a[19].order[2]=a[19].order[3]=0;
	a[19].norder=2;
	a[19].nbond=7;
	a[19].rb=1.05;
	a[19].type=4;
	a[19].bd[0]=2;
	a[19].ang0=104.5*f;
	a[19].atm="S";

	a[20].id=20;
	a[20].name="S(=)";
	a[20].order[0]=2;
	a[20].order[1]=a[20].order[2]=a[20].order[3]=0;
	a[20].nbond=4;
	a[20].norder=1;
	a[20].rb=1.05;
	a[20].type=7;
	a[20].bd[1]=1;
	a[20].ang0=0;
	a[20].atm="S";

	a[21].id=21;
	a[21].name="P(-)(-)(-)";
	a[21].order[0]=a[21].order[1]=a[21].order[2]=1;
	a[21].order[3]=0;
	a[21].nbond=10;
	a[21].norder=3;
	a[21].rb=1.07;
	a[21].type=2;
	a[21].bd[0]=3;
	a[21].ang0=107.5*f;
	a[21].atm="P";

	a[22].id=22;
	a[22].name="P(=)(-)";
	a[22].order[0]=2;
	a[22].order[1]=1;
	a[22].order[2]=a[22].order[3]=0;
	a[22].nbond=7;
	a[22].norder=2;
	a[22].rb=1.07;
	a[22].type=4;
	a[22].bd[0]=1;
	a[22].bd[1]=1;
	a[22].ang0=115.0*f;
	a[22].atm="P";
	
	a[23].id=23;
	a[23].name="P(#)";
	a[23].order[0]=3;
	a[23].order[1]=a[23].order[2]=a[23].order[3]=0;
	a[23].nbond=4;
	a[23].norder=1;
	a[23].rb=1.07;
	a[23].type=7;
	a[23].bd[2]=1;
	a[23].ang0=0;
	a[23].atm="P";

	a[24].id=24;
	a[24].name="[F-]";
	for (i=0;i<4;i++) a[24].order[i]=0;
	a[24].nbond=4;
	a[24].norder=0;
	a[24].rb=0.72;
	a[24].type=0;
	a[24].chg=-1;
	a[24].atm="F";

	a[25].id=25;
	a[25].name="[Cl-]";
	for (i=0;i<4;i++) a[25].order[i]=0;
	a[25].nbond=5;
	a[25].norder=0;
	a[25].rb=0.99;
	a[25].type=0;
	a[25].chg=-1;
	a[25].atm="Cl";

	a[26].id=26;
	a[26].name="[Br-]";
	for (i=0;i<4;i++) a[26].order[i]=0;
	a[26].nbond=5;
	a[26].norder=0;
	a[26].rb=1.20;
	a[26].type=0;
	a[26].chg=-1;
	a[26].atm="Br";

	a[27].id=27;
	a[27].name="[I-]";
	for (i=0;i<4;i++) a[27].order[i]=0;
	a[27].nbond=4;
	a[27].norder=0;
	a[27].rb=1.39;
	a[27].type=0;
	a[27].chg=-1;
	a[27].atm="I";

	a[28].id=28;
	a[28].name="[OH-]";
	for (i=0;i<4;i++) a[28].order[i]=0;
	a[28].nbond=5;
	a[28].norder=0;
	a[28].chg=-1;
	a[28].nh=1;
	a[28].rb=0.74;
	a[28].type=0;	
	a[28].atm="O";

	a[29].id=29;
	a[29].name="[OH0-](-)";
	a[29].order[0]=1;
	a[29].nbond=9;
	a[29].norder=1;
	a[29].chg=-1;
	a[29].rb=0.74;
	a[29].type=7;
	a[29].bd[0]=1;
	a[29].index=7;	
	a[29].atm="O";

	a[30].id=30;
	a[30].name="P(-)(-)(-)(-)(-)(-)";
	for (i=0;i<6;i++) a[30].order[i]=1;
	a[30].nbond=19;
	a[30].norder=6;
	a[30].chg=-1;
	a[30].rb=1.07;
	a[30].type=6;
	a[30].bd[0]=6;
	a[30].ang0=1.82557;
	a[30].atm="P";

	a[31].id=31;
	a[31].name="P(-)(-)(-)(-)(-)";
	for (i=0;i<5;i++) a[31].order[i]=1;
	a[31].nbond=16;
	a[31].norder=5;
	a[31].chg=0;
	a[31].rb=1.07;
	a[31].type=5;
	a[31].bd[0]=5;
	a[31].ang0=90.0*f;
	a[31].atm="P";

	a[32].id=32;
	a[32].name="P(=)(-)(-)(-)";
	a[32].order[0]=2;
	for (i=1;i<4;i++) a[32].order[i]=1;
	a[32].nbond=13;
	a[32].norder=4;
	a[32].chg=0;
	a[32].rb=1.07;
	a[32].type=1;
	a[32].bd[0]=3;
	a[32].bd[1]=1;
	a[32].ang0=109.5*f;
	a[32].atm="P";

	a[33].id=33;
	a[33].name="S(=)(=)(-)(-)";
	a[33].order[0]=a[33].order[1]=2;
	a[33].order[2]=a[33].order[3]=1;
	a[33].nbond=13;
	a[33].norder=4;
	a[33].chg=0;
	a[33].rb=1.05;	
	a[33].type=1;
	a[33].bd[0]=2;
	a[33].bd[1]=2;
	a[33].ang0=109.5*f;
	a[33].atm="S";

	a[34].id=34;
	a[34].name="S(=)(-)(-)";
	a[34].order[0]=2;
	a[34].order[1]=a[34].order[2]=1;
	a[34].norder=3;
	a[34].nbond=10;
	a[34].chg=0;
	a[34].rb=1.05;
	a[34].type=3;
	a[34].bd[0]=2;
	a[34].bd[1]=1;
	a[34].ang0=120.0*f;
	a[34].atm="S";

	a[35].id=35;
	a[35].name="[NH0-](-)(-)";
	a[35].order[0]=a[35].order[1]=1;
	a[35].nbond=12;
	a[35].norder=2;
	a[35].chg=-1;
	a[35].rb=0.74;
	a[35].type=4;
	a[35].bd[0]=2;
	a[35].chg=-1;
	a[35].index=7;
	a[35].ang0=104.5*f;
	a[35].atm="N";	
	return;
}
