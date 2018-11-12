#include "syspara.h"

void function(double x[],double f[],double t)
{

	int i;
	comp_reversal_potential(x);
	comp_CaMK(x);
	comp_ina(x);
	comp_ito(x);
	comp_ical(x);
	comp_ikr(x);
	comp_iks(x);
	comp_ik1(x);
	comp_inaca(x);
	comp_inak(x);
	comp_ipca(x);
	comp_ikb(x);
	comp_icab(x);
	comp_inab(x);
	//comp_CaMK(x);
	comp_diff(x);
	comp_jrel(x);
	comp_jup(x);
	comp_jtr(x);
	comp_concentration(x);
	
	var.Ina_i_total = ina.total + inab.na + 3.0*inak.inak + 3.0*var.inaca_i;
	var.Ina_ss_total = ical.icana + 3.0*var.inaca_ss;
	var.Ik_i_total = ito.ik + ikr.ik + iks.ik + ik1.ik + ikb.k - 2.0*inak.inak + var.Istim;
	var.Ik_ss_total = ical.icak;
	var.Ica_i_total = ipca.ca + icab.ca - 2.0*var.inaca_i;
	var.Ica_ss_total = ical.ica - 2.0*var.inaca_ss;
	var.Itotal = ina.total + ito.ik + ical.ica + ical.icana + ical.icak + ikr.ik + iks.ik + ik1.ik 
					+ var.inaca + inak.inak + inab.na + icab.ca + ikb.k + ipca.ca + var.Istim;

	f[0] = -var.Itotal;
	//Fast sodium current
	f[1] = (ina.mss - x[1])/ina.taum; // m
	f[2] = (ina.hss - x[2])/ina.tauh; // h_fast
	f[3] = (ina.jss - x[3])/ina.tauj; // j
	f[4] = (ina.hCaMKss - x[4])/ina.tauh_CaMK_slow; // h_CaMK_slow
	f[5] = (ina.jss - x[5])/ina.tauj_CaMK; // j_CaMK
	//late sodium current
	f[6] = (ina.mlss - x[6])/ina.tauml; // ml
	f[7] = (ina.hlss - x[7])/ina.tauhl; // hl
	f[8] = (ina.hlCaMKss - x[8])/ina.tauhl_CaMK; // hl
	//Transient outward current
	f[9] = (ito.ass - x[9])/ito.taua;
	f[10] = (ito.iss - x[10])/ito.taui_fast;
	f[11] = (ito.iss - x[11])/ito.taui_slow;
	f[12] = (ito.aCaMKss - x[12])/ito.taua;
	f[13] = (ito.iss - x[13])/ito.taui_CaMK_fast;
	f[14] = (ito.iss - x[14])/ito.taui_CaMK_slow;
	// LTCC
	f[15] = (ical.dss - x[15])/ical.taud;
	f[16] = (ical.fss - x[16])/ical.tauf_fast;
	f[17] = (ical.fss - x[17])/ical.tauf_slow;
	f[18] = (ical.fss - x[18])/ical.taufca_fast;
	f[19] = (ical.fss - x[19])/ical.taufca_slow;
	f[20] = (ical.fss - x[20])/ical.taujca;
	f[21] = (ical.fss - x[21])/ical.tauf_CaMK_fast;
	f[22] = (ical.fss - x[22])/ical.taufca_CaMK_fast;
	f[23] = ical.alpha_n*ical.kp2n - x[23]*ical.km2n;
	// Ikr
	f[24] = (ikr.xrss - x[24])/ikr.tauxr_fast;
	f[25] = (ikr.xrss - x[25])/ikr.tauxr_slow;
	// Iks
	f[26] = (iks.xs1ss - x[26])/iks.tauxs1;
	f[27] = (iks.xs1ss - x[27])/iks.tauxs2;
	// Ik1
	f[28] = (ik1.k1ss - x[28])/ik1.tauk1;
	// CaMK
	f[29] = CaMK.a*CaMK.bound*(CaMK.bound+x[29]) - CaMK.b*x[29];
	// Jrel
	f[30] = (jrel.NPss - x[30])/jrel.tau_NP; 
	f[31] = (jrel.CaMKss - x[31])/jrel.tau_CaMK; 
	// [Na]i
	f[32] = -var.Ina_i_total*var.vr1 + jdiff.na*var.vr2;
	// [Na]ss
	f[33] = -var.Ina_ss_total*var.vr3 - jdiff.na;
	// [K]i
	f[34] = -var.Ik_i_total*var.vr1 + jdiff.k*var.vr2;
	// [K]ss
	f[35] = -var.Ik_ss_total*var.vr3 - jdiff.k;
	// [Ca]i
	f[36] = var.b_Ca_i*(-var.Ica_i_total*var.vr4 - jup.ca*var.vr5 + jdiff.ca*var.vr2);
	// [Ca]ss
	f[37] = var.b_Ca_ss*(-var.Ica_ss_total*var.vr6 +jrel.ca*var.vr7 - jdiff.ca);
	// [Ca]nsr
	f[38] = jup.ca - jtr.ca*var.vr8;
	// [Ca]jsr
	f[39] = var.b_Ca_jsr*(jtr.ca-jrel.ca);

	//printf("NPss=%lf,NPp=%lf\n",jrel.NPss,jrel.CaMKss);
	//for(i=0;i<NN;i++){
	//	printf("x[%d]=%e\n",i,f[i]);
	//}
}
