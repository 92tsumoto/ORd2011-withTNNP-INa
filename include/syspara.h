//#ifndef __SYSPARA_H_INCLUDE 
//#define __SYSPARA_H_INCLUDE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "/home/tsumoto/lib/xhplot.h"

#define NN 24
#define BUF 100
#define NUM 20

#define R 8314.472
#define F 96485.33771638995
#define T 310

struct varstruct {

    int datas;
    int line_wid[NUM];
	
	int n;
    double vm,vs,ks,vd,k1,k2,tmp,Istim;
	double dvm,dvs,dks,dvd,dk1,dk2,dtmp,dIstim;

	// Cell Geometry
	double length,a,vcell,ageo,acap,vmyo,vmito,vsr,vnsr,vjsr,vcleft;

	// Na Ion concentration
	//double nai, nao;
	// Ca Ion concentration
	//double cai, cao;
	// K Ion concentration
	//double ki, ko;

	// Ion Valences 
	double zna,zk,zca;
	
	// Fast sodium current
	double gna,Ena,am,bm,ah,bh,aj,bj,ina;

	// Change rates for K channel conductances
	double ikrf,iksf,ikxf,itof;

	// L-type calcium current
	double ilcatot,ibarca,ibarna,ibark;
	double ilcarat;
	double ilca,ilcana,ilcak;
	double dss,taud;                  // Steady-state value of activation gate d and time const
	double fss1,fss2,tauf1,tauf2;     // Steady-state value of inactivation gate f and time const
	double fcass1,fcass2,fcatau1,fcatau2; // Ca dependant inactivation gate
	double kmca,pca,gacai,gacao;
	double pna,ganai,ganao;           // Permiability of membrane to Na (cm/s)
	double pk,gaki,gako;              // Permiability of membrane to K (cm/s)
	double ratgca,fcarat,aCDI,bCDI;

	// T-type calcium current
	double gcat,Eca,bss,taub,gss,taug,icat;

	// Rapidly activating potassium current
	double gkr,Ekr,xrss,tauxr,r,ikr,gkr_max;
	
	// Slowly activating potassium current
	double gks,Eks,xs1ss,tauxs1,xs2ss,tauxs2;
	double prnak,iks,gks_max;

	// Potassium current (time-dependent)
	double gki,Eki,aki,bki,kin,iki;

	// Plateau Potassium current
	double gkp,Ekp,kp,ikp,gkp_max;

	// Ultra rapid Potassium current
	double gkur,ikur;

	// Na-activated potassium channel
	double gkna,pona,pov,Ekna,kdkna,ikna;

	// ATP-sensitive potassium channel
	double gkatp,gkbaratp,patp,natp,Ekatp;
	double nicholsarea,atpi,hatp,katp,ikatp;

	// Ito Transient Outward Current
	// (Dumaine et al. Circ Res 1999;85:803-809)
	double gitodv,Ekdv,rvdv,azdv,bzdv;
	double aydv,bydv,ito,gitodv_max;

	 // Sodium-Calcium Exchanger V-S
	double c1,c2,gammas,inaca;

	// Sodium-Potassium Pump
	double fnak,sigma,ibarnak,kmnai,kmko,inak;
	
	// Nonspecific Ca-activated Current 
	double ibarnsna,ibarnsk,pnsca,kmnsca,insna,insk;

	// Sarcolemmal Ca Pump
	double ibarpca,kmpca,ipca;

	// Ca Background Current
	double gcab,Ecan,icab;

	// Na Background Current
	double gnab,Enan,inab;

	// Total Ion currents 
	double INa_total, IK_total, ICa_total;
	
	// difference total Ion current 
	double dICa_total;
	
	// NSR Ca Ion Concentration Changes
	double iup,ileak;      // Ca uptake from myo. to NSR (mM/ms)
	double kmup,iupbar,nsrbar;

	// JSR Ca Ion Concentration Changes
	double Irel_cicr,Irel_jsr_overload;
	double dICa_total_new,dCa_ion,ICa_total_old,boolien,gmaxrel;
	double jsr_new,grelbarjsrol,tauon,tauoff;
	double csqnbar,kmcsqn,t_cicr,t_overload;
	double bjsr,cjsr,csqn,csqnth,swspontan;
	double djsr,jsr;

	// test variable
	double dt;

	// Translocation of Ca Ions from NSR to JSR
	double tautr,itr;

	// Ca concentration
	double cmdnbar,trpnbar,kmtrpn,kmcmdn;
	double trpn,cmdn;
	double Ca_total,gpig,bmyo,cmyo,dmyo;

	// Extracellular ion concentrations
	double Na_bulk,K_bulk,Ca_bulk;
	double tau_diff;

	// Base Currnt Stimulus
	double Istim_base;

	// Sttimulus parameters
	double BCL;  // Base cycle length = stimulus period
	int beat; // Number of stimulus


	// debug variable
	double ca_pre,dca_now;

    int m;
    int l;

    double x0[NUM][NN];
    double tsign[NUM];
    double tend[NUM];

    int pflag;
    int write, graph;
    int write0;
    int half;

} var;


void eular(int n,double h,double x[],double t);

void function(double x[],double f[],double t);

void input_para(FILE *);

void eventloop(FILE *, int *mode, int *P, double m[]);

void orbit(int *mode, double m[]);

void draw_p(int *mode, int P, double x[]);

void mouse(int *mode, double x[]);

void comp_ina(double x[]);
void comp_icat(double x[]);
void comp_ical(double x[]);
void comp_ikr(double x[]);
void comp_iki(double x[]);
void comp_iks(double x[]);
void comp_ikp(double x[]);
void comp_ikur(double x[]);
void comp_ikna(double x[]);
void comp_ikatp(double x[]);
void comp_ito(double x[]);
void comp_inaca(double x[]);
void comp_inak(double x[]);
void comp_insca(double x[]);
void comp_ipca(double x[]);
void comp_icab(double x[]);
void comp_inab(double x[]);
void comp_iconcent (double x[]);
void comp_iconcent2 (double x[]);
void conc_nsr(double x[]);
void conc_jsr(double x[]);
void conc_itr (double x[]);
void conc_cai (double x[]);
void conc_cleft (double x[]);

void main(int argc,char **argv);
