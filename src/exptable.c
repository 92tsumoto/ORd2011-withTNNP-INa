#include "syspara.h"

void make_ExpTable()
{

	int vindex,kiindex;
	double v,ki;
	double am,bm,ah,bh,aj,bj;
    
	for(vindex=0;vindex<VNMAX;vindex++){

        v = (double)vindex/dvm-200.0;
		
        /** for original ina **/
		//ina.Tmss[vindex] = 1.0/(1.0+exp(-(v+39.57)/9.871));
		//ina.Ttaum[vindex] = 1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
        
		/** for ina TNNP **/
		ina.Tmss[vindex] = 1.0/((1.0+exp(-(v+56.86)/9.03))*(1.0+exp(-(v+56.86)/9.03)));
		am=1.0/(1.+exp((-60.0-v)/5.0));
		bm=0.1/(1.0+exp((v+35.0)/5.0))+0.10/(1.0+exp((v-50.0)/200.0));
		ina.Ttaum[vindex] = am*bm;

        /** for original ina **/
		//ina.Thss[vindex] = 1.0/(1.0+exp((v+82.90)/6.086));
		//ina.Ttauh_fast[vindex] = 1.0/(1.432E-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));
		//ina.Ttauh_slow[vindex] = 1.0/(0.009764*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
		
		/** for ina TNNP **/
		ina.Thss[vindex] = 1.0/((1.0+exp((v+71.55)/7.43))*(1.0+exp((v+71.55)/7.43)));
		if(v>=-40.0){
			ah = 0.0;
			bh = 0.77/(0.13*(1.0+exp(-(v+10.66)/11.1)));
		} else {
			ah = 0.057*exp(-(v+80.0)/6.8);
			bh = 2.7*exp(0.079*v)+3.1E+5*exp(0.3485*v);
		}
		ina.Ttauh[vindex] = 1.0/(ah+bh);

        /** for original ina **/
		//ina.Ttauj[vindex] = 2.038+1.0/(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45));

		/** for ina TNNP **/
		ina.Tjss[vindex] = 1.0/((1.0+exp((v+71.55)/7.43))*(1.0+exp((v+71.55)/7.43)));
		if(v>=-40.0){
			aj = 0.0;
			bj = 0.6*exp(0.057*v)/(1.0+exp(-0.1*(v+32.0)));
		} else {
			aj = (-2.5428E+4*exp(0.2444*v)-6.948E-6*exp(-0.04391*v))*(v+37.78)/(1.0+exp(0.311*(v+79.23)));
			bj = 0.02424*exp(-0.01052*v)/(1.0+exp(-0.1378*(v+40.14)));
		}
		ina.Ttauj[vindex] = 1.0/(aj+bj);


		ina.ThCaMKss[vindex] = 1.0/(1.0+exp((v+89.1)/6.086));

		// for inal 
		ina.Tmlss[vindex] = 1.0/(1.0 + exp(-(v+42.85)/5.264));
		ina.Ttauml[vindex] = 1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
		ina.Thlss[vindex] = 1.0/(1.0 + exp((v + 87.61)/7.488));
		ina.ThlCaMKss[vindex] = 1.0/(1.0 + exp((v+93.81)/7.488));

		// ito
		ito.Tass[vindex] = 1.0/(1.0+exp(-(v-14.34)/14.82));
		ito.Ttaua[vindex] = 1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));

		ito.Tiss[vindex] = 1.0/(1.0+exp((v+43.94)/5.711));
		ito.Ttaui_fast[vindex] = 4.562+1.0/(0.3933*exp(-(v+100.0)/100.0)+0.08004*exp((v+50.0)/16.59));
		ito.Ttaui_slow[vindex] = 23.62+1.0/(0.001416*exp(-(v+96.52)/59.05)+1.7808E-8*exp((v+114.1)/8.079));
		ito.TAi_fast[vindex] = 1.0/(1.0+exp((v-213.6)/151.2));
		ito.TaCaMKss[vindex] = 1.0/(1.0+exp(-(v-24.34)/14.82));
		ito.TdeltaCaMK_dev[vindex] = 1.354+1.0E-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
		ito.TdeltaCaMK_rec[vindex] = 1.0-0.5/(1.0+exp((v+70.0)/20.0));
		if(var.celltype==3){
			ito.Tdepi[vindex] = 1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
		} 

		// for ical
		ical.Tdss[vindex] = 1.0/(1.0+exp(-(v+3.940)/4.230));
		ical.Ttaud[vindex] = 0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));

		ical.Tfss[vindex] = 1.0/(1.0+exp((v+19.58)/3.696));
		ical.Ttauf_fast[vindex] = 7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
		ical.Ttauf_slow[vindex] = 1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));

		ical.Ttaufca_fast[vindex] = 7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
		ical.Ttaufca_slow[vindex] = 100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
		ical.TAfca_fast[vindex] = 0.3+0.6/(1.0+exp((v-10.0)/10.0));

		ical.Texp_Ca[vindex] = exp(v/var.RTon2F);
		ical.Texp_Na[vindex] = exp(v/var.RTonF);
		ical.Texp_K[vindex] = exp(v/var.RTonF);
		
		// for ikr 
		ikr.Txrss[vindex] = 1.0/(1.0+exp(-(v+8.337)/6.789));
		ikr.Ttauxr_fast[vindex] = 12.98+1.0/(0.3652*exp((v-31.66)/3.869)+4.123E-5*exp(-(v-47.78)/20.38));
		ikr.Ttauxr_slow[vindex] = 1.865+1.0/(0.06629*exp((v-34.70)/7.355)+1.128E-5*exp(-(v-29.74)/25.94));
		ikr.TAxr_fast[vindex] = 1.0/(1.0+exp((v+54.81)/38.21));
		ikr.Trkr[vindex] = 1.0/(1.0+exp((v+55.0)/75.0))*1.0/(1.0+exp((v-10.0)/30.0));

		// for iks 

		iks.Txs1ss[vindex] = 1.0/(1.0+exp(-(v+11.60)/8.932));
		iks.Ttauxs1[vindex] = 817.3+1.0/(2.326E-4*exp((v+48.28)/17.80)+0.001292*exp(-(v+210.0)/230.0));
		iks.Ttauxs2[vindex] = 1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp(-(v+66.54)/31.0));

		// ik1 
		ik1.Tk1ss[vindex] = 1.0/(1.0+exp(-(v+2.5538*var.ko+144.59)/(1.5692*var.ko+3.8115)));
		ik1.Ttauk1[vindex] = 122.2/(exp(-(v+127.2)/20.36)+exp((v+236.8)/69.33));
		ik1.Trk1[vindex] = 1.0/(1.0+exp((v+105.8-2.6*var.ko)/9.493));

		// inaca
		var.Thca[vindex] = exp(var.qca*v/var.RTonF);
		var.Thna[vindex] = exp(var.qna*v/var.RTonF);

		// inak 
		inak.Tknai[vindex] = inak.ko_nai*exp((inak.delta*v*F)/(3.0*R*T));
		//inak.Tknao[vindex] = inak.ko_nao*exp(((1.0-inak.delta)*v*F)/(3.0*R*T));
		inak.Tknao[vindex] = inak.ko_nao*exp(((1.0+inak.delta)*v*F)/(3.0*R*T));

		// ikb
		ikb.Txkb[vindex] = 1.0/(1.0+exp(-(v-14.48)/18.34));

		// icab
		icab.Texp[vindex] = exp(v/var.RTon2F);

		// inab
		inab.Texp[vindex] = exp(v/var.RTonF);
	}

}
