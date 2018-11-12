//#ifndef __SYSPARA_H_INCLUDE 
//#define __SYSPARA_H_INCLUDE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mkl.h"
#include "./include/xhplot.h"

#define NN 40
#define BUF 100
#define NUM 60

//#define R 8314.472		// J/kmol/K
//#define F 96485.33771638995	// C/mol
//#define T 310.0		// K
#define R 8314.0		// J/kmol/K
#define F 96485.0	// C/mol
#define T 310.0		// K

#define VNMAX 400*5+1
#define dvm 5

struct varstruct {

    int datas;
    int line_wid[NUM];
	
	int n;
    double Istim,dIstim;
	double rategna,rategca;

	// An invariant constant
	double RTonF,RTon2F;

	// Cell tupe
	int celltype;

	// Cell Geometry
	double length,a,vcell,ageo,acap;
	double vmyo,vmito,vsr,vnsr,vjsr,vcleft,vss;
	double vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8;

	// Ion Valences 
	double zna,zk,zca;

	// Reversal potential
	double Ena,Ek,Eks;
	double prnak;
			
	// Sodium-Calcium Exchanger V-S
	double hca,hna;
	double *Thca,*Thna;
	double kna1,kna2,kna3,kasym;
	double omega_na,omega_ca,omega_naca;
	double kca_on,kca_off,qna,qca;
	double km_ca_act;
	double inaca_i,inaca_ss,inaca;
	double Gnaca;

	// Total Ion currents 
	double Ina_i_total, Ina_ss_total;
	double Ik_i_total, Ik_ss_total;
	double Ica_i_total,Ica_ss_total;
	double Itotal;

	// Ion concentration and Buffers
	double cmdnbar,kmcmdn;
	double trpnbar,kmtrpn;
	double bsrbar,kmbsr;
	double bslbar,kmbsl;
	double csqnbar,kmcsqn;
	double b_Ca_i,b_Ca_ss,b_Ca_jsr;

	// Extracellular ion concentrations
	double nao,ko,cao;

	// Base Currnt Stimulus
	double Istim_base;

	// test variable
	double dt;
	// Sttimulus parameters
	double BCL;  // Base cycle length = stimulus period
	int beat; // Number of stimulus

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

// Fast and Late sodium currnets
struct inastruct {

	double *Tmss,*Ttaum,*Thss,*Ttauh,*Ttauh_fast,*Ttauh_slow,*Tjss,*Ttauj;
	double *ThCaMKss;
	double Gna_fast,fast;
	double mss,taum,hss,tauh,tauh_fast,tauh_slow,jss,tauj;
	double Ah_fast,Ah_slow;
	double hCaMKss,tauh_CaMK_slow,jCaMKss,tauj_CaMK;
	double Ah_CaMK_fast,Ah_CaMK_slow;
	double h,hCaMK;
	double fast_pCaMK;

	// Late sodium currnet
	double *Tmlss,*Ttauml,*Thlss,*ThlCaMKss;
	double Gna_late,late;
	double mlss,tauml,hlss,tauhl,hlCaMKss,tauhl_CaMK;
	double late_pCaMK;

	// ina
	double total;
} ina;

// Transient Outward Current (Ito)
struct itostruct {

	double ik,Gto;
	double ass,taua,iss,taui_fast,taui_slow,Ai_fast,Ai_slow,i,depi;
	double *Tass,*Ttaua,*Tiss,*Ttaui_fast,*Ttaui_slow,*TAi_fast,*Tdepi;
	double aCaMKss,taua_CaMK,iCaMKss,deltaCaMK_dev,deltaCaMK_rec;
	double *TaCaMKss,*TdeltaCaMK_dev,*TdeltaCaMK_rec;
	double taui_CaMK_fast,taui_CaMK_slow,Ai_CaMK_fast,Ai_CaMK_slow,iCaMK;
	double pCaMK;

} ito;

// L-type Calcium channel current (IcaL)
struct icalstruct {

	double dss,taud,fss,tauf_fast,tauf_slow;
	double *Tdss,*Ttaud,*Tfss,*Ttauf_fast,*Ttauf_slow;
	double f,Af_fast,Af_slow;

	double fcass,taufca_fast,taufca_slow,Afca_fast,Afca_slow,fca; 
	double *Ttaufca_fast,*Ttaufca_slow,*TAfca_fast; 

	double jcass,taujca; 

	double f_CaMKss,f_CaMK,tauf_CaMK_fast,Af_CaMK_fast,Af_CaMK_slow,f_CaMK_slow; 

	double fca_CaMKss,fca_CaMK,taufca_CaMK_fast,Afca_CaMK_fast,Afca_CaMK_slow,fca_CaMK_slow; 

	double kmn,kp2n,km2n,alpha_n;
	double pca,gacai,gacao,phi_ca,ibarcal;
	double exp_Ca,exp_Na,exp_K;
	double *Texp_Ca,*Texp_Na,*Texp_K;
	double pcana,ganai,ganao,phi_na,ibarcana; 
	double pcak,gaki,gako,phi_k,ibarcak;

	double pca_CaMK,pcana_CaMK,pcak_CaMK,ibarcal_CaMK,ibarcana_CaMK,ibarcak_CaMK;
	double pCaMK;
	double ica,icana,icak;
	double tmp;
} ical;

// Rapid activating potassium current (Ikr)
struct ikrstruct {

	double ik,Gkr,rategkr;
	double xr,xrss,tauxr_fast,tauxr_slow,Axr_fast,Axr_slow,rkr;
	double *Txrss,*Ttauxr_fast,*Ttauxr_slow,*TAxr_fast;
	double *Trkr;

} ikr;

// Slowlactivating potassium current (Iks)
struct iksstruct {

	double ik,Gks,KsCa;
	double xs1ss,tauxs1,xs2ss,tauxs2;
	double *Txs1ss,*Ttauxs1,*Ttauxs2;
		
} iks;

// Inward rectifier potassium current (Ik1)
struct ik1struct {

	double ik,Gk1,rategk1;
	double k1ss,tauk1,rk1;
	double *Tk1ss,*Ttauk1,*Trk1;

} ik1;

struct ncxistruct {

	double h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12;
	double k1,k2,k3,k31,k32,k4,k41,k42;
	double k5,k6,k7,k8;
	double x1,x2,x3,x4;
	double E1,E2,E3,E4;
	double allo;
	double jnaca_na,jnaca_ca;

} ncxi;

struct ncxssstruct {

	double h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12;
	double k1,k2,k3,k31,k32,k4,k41,k42;
	double k5,k6,k7,k8;
	double x1,x2,x3,x4;
	double E1,E2,E3,E4;
	double allo;
	double jnaca_na,jnaca_ca;

} ncxss;

// Sodium-Potassium Pump
struct inakstruct {

	double kp1,km1,kp2,km2;
	double kp3,km3,kp4,km4;
	double ko_nai,ko_nao,delta;
	double knai,knao,*Tknai,*Tknao;
	double kki,kko,MgADP,MgATP,k_MgATP,H,SigP;
	double HP,nap,kp,P;
	double a1,a2,a3,a4;
	double b1,b2,b3,b4;
	double x1,x2,x3,x4;
	double E1,E2,E3,E4;
	double G,jna,jk,inak;

} inak;

// Sarcolemmal Ca Pump
struct ipcastruct {

	double G,km,ca;

} ipca;

// Na Background Current
struct inabstruct {

	double pnab,na;
	double exp,*Texp;


} inab;

// K Background Current
struct ikbstruct {

	double G,k;
	double xkb,*Txkb;

} ikb;

// Ca Background Current
struct icabstruct {

	double pcab,gacai,gacao,ca;
	double exp,*Texp;

} icab;

// Ca/Calmodulin dependent protein kinese (CaMK)
struct CaMKstruct{
	
	double a,b,z,bound,active;
	double Km;

} CaMK;

// Ca/Calmodulin dependent protein kinese (CaMK)
struct CaMstruct{
	
	double Km;

} CaM;

// SR calcium release flux, via RyR (Jrel)
struct jrelstruct {

	double b_tau,a,b_tau_CaMK,a_CaMK;
	double NPss,tau_NP;
	double CaMKss,tau_CaMK;
	double pCaMK,ca;
	double p;

} jrel;

// Calcium uptake via SERCA pump
struct jupstruct {

	double np;
	double dKm_PLB,dCaMK,CaMK;
	double pCaMK,ca,leak;
	double p;

} jup;

// diffusion flux
struct jdiffstruct {

	double tau_na,tau_k,tau_ca;
	double na,k,ca;

} jdiff;

// Translocation of Ca Ions from NSR to JSR
struct jtrstruct {

	double tau,ca;

} jtr;

void val_consts(FILE *);
void make_ExPTable();

void eular(int n,double h,double x[],double t);
void function(double x[],double f[],double t);
void input_para(FILE *);
//void initial_mem(int tMAX);
void initial_mem();
void closed_mem();

void eventloop(FILE *, int *mode, int *P, double m[]);
void orbit(int *mode, double m[], double x2);
void draw_p(int *mode, int P, double x[], double x2);
void mouse(int *mode, double x[], double x2);

void comp_reversal_potential(double x[]);
void comp_ina(double x[]);
//void comp_inal(double x[]);
void comp_ito(double x[]);
void comp_ical(double x[]);
void comp_ikr(double x[]);
void comp_iks(double x[]);
void comp_ik1(double x[]);
void comp_inaca(double x[]);
void comp_inak(double x[]);
void comp_ipca(double x[]);
void comp_ikb(double x[]);
void comp_icab(double x[]);
void comp_inab(double x[]);
void comp_CaMK(double x[]);
void comp_diff(double x[]);
void comp_jrel(double x[]);
void comp_jup(double x[]);
void comp_jtr (double x[]);
void comp_concentration (double x[]);
//void comp_iconcent (double x[]);
//void comp_iconcent2 (double x[]);
//void conc_nsr(double x[]);
//void conc_jsr(double x[]);
//void conc_itr (double x[]);
//void conc_cai (double x[]);
//void conc_cleft (double x[]);

main(int argc,char **argv);
