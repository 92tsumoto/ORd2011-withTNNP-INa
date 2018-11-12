
#include "syspara.h"

void val_consts(FILE *fp1)
{
	int i,w;
	double v_old,dvdt,dvdt_new;

	// Cell Geometry */
	// Forbes M, Sperelakis N. Ultrastructure of mammalian cardiac muscle. 
	// In: Sperelakis N, ed. Physiology and pathophysiology of the heart.
	// 2nd edition. Boston, MA: Kluwer Academic; 1989:3-41.
		//var.length = 0.01;       // Length of the cell (cm)
		//var.a = 0.0011;     // Radius of the cell (cm)
		var.length = 0.0135;       // Length of the cell (cm)
		var.a = 0.001;     // Radius of the cell (cm)
		var.vcell = 1000*M_PI*var.a*var.a*var.length; // Cell Volume:3.801e-5 (uL)
		var.ageo = 2*M_PI*var.a*var.a+2.0*M_PI*var.a*var.length;  // eometric membrane area: 7.671e-5 (cm^2)
		var.acap = var.ageo*2;          // Capacitive membrane area: 1.534e-4 cm^2 (cm^2)
		var.vmyo = var.vcell*0.68;      // Myoplasm volume (uL) = 68% for Cell volume
		var.vnsr = var.vcell*0.0552;    // NSR volume (uL)
		var.vjsr = var.vcell*0.0048;    // JSR volume (uL)
		var.vss  = var.vcell*0.02;    	// Subspace compartment (representing submembrane space near t-tubles)
		//var.vcleft = var.vcell*0.12/0.88;  // Cleft volume (uL)
		//var.vmito = var.vcell*0.26;     // Mitochondria volume (uL) = 26% for cell volume
		//var.vsr = var.vcell*0.06;       // SR volume (uL)
		var.vr1 = var.acap/(F*var.vmyo);
		var.vr2 = var.vss/var.vmyo;
		var.vr3 = var.acap/(F*var.vss);
		var.vr4 = var.acap/(2.0*F*var.vmyo);
		var.vr5 = var.vnsr/var.vmyo;
		var.vr6 = var.acap/(2.0*F*var.vss);
		var.vr7 = var.vjsr/var.vss;
		var.vr8 = var.vjsr/var.vnsr;


	// Ion Valences
		var.zna = 1.0;  // Na valence
		var.zk = 1.0;   // K valence
		var.zca = 2.0;  // Ca valence

	// invariant constant
		var.RTonF = R*T/F;
		var.RTon2F = R*T/(var.zca*F);

	// Extracellular Concentrations
		var.nao = 140.0;     // Initial Bulk Medium Na (mM)
		var.ko = 5.4;      // Initial Bulk Medium K (mM)
		var.cao = 1.8;     // Initial Bulk Medium Ca (mM)
		var.prnak = 0.01833;     // Na/K Permiability Ratio

	// Fast sodium current
		//ina.Gna_fast = 75.0;	// (mS/uF).
		//ina.Gna_fast = 14.868;	// (mS/uF).
		ina.Gna_fast = var.rategna*11.0;	// (mS/uF).
		//ina.Gna_fast = 0.33333333*11.0;	// (mS/uF).
		ina.Ah_fast = 0.99;
		ina.Ah_slow = 0.01; 	// (1.0-ina.Ah_fast)
		//ina.Ah_CaMK_fast = ina.Ah_fast;
		//ina.Ah_CaMK_slow = ina.Ah_slow;
		
	// Late sodium current
		ina.tauhl = 200.0;				// (ms).
		ina.tauhl_CaMK = 3.0*ina.tauhl;	// 600 (ms).
		if(var.celltype==1){	// Endo.
			ina.Gna_late = 0.0075;		
		} else if(var.celltype==2){ // Mid.
			ina.Gna_late = 0.0075;	
		} else if(var.celltype==3){
			ina.Gna_late = 0.0075*0.6;	// Epi.
		}
		
	// Transient outward current
		if(var.celltype==1){ 	// Endo.
			ito.Gto = 0.02;	// (mS/uF).
		} else if(var.celltype==2){	// Mid.
			ito.Gto = 4.0*0.02;
		} else if(var.celltype==3){	// Epi.
			ito.Gto = 4.0*0.02;
		}

	// L-type calcium current
		ical.Af_fast = 0.6;
		ical.Af_slow = 1.0 - ical.Af_fast;
		ical.taujca = 75.0; // (ms).
		//ical.Af_CaMK_fast = ical.Af_fast;
		//ical.Af_CaMK_slow = ical.Af_slow;

		ical.kmn = 0.002;
		ical.kp2n = 1000.0;

		ical.gacai = 1.0;         					// Activity coefficient of Ca
		ical.gacao = 0.341;     					// Activity coefficient of Ca
		ical.ganai = 0.75;      					// Activity coefficient of Na
		ical.ganao = 0.75;      					// Activity coefficient of Na
		ical.gaki = 0.75;       					// Activity coefficient of K
		ical.gako = 0.75;       					// Activity coefficient of K
	
		if(var.celltype == 1){	// Endo
			ical.pca = var.rategca*0.0001;      					// Permiability of membrane to Ca (cm/s)
			ical.pcana = 0.00125*ical.pca;   			// Permiability of membrane to Na (cm/s)
			ical.pcak = 3.574E-4*ical.pca;  			// Permiability of membrane to K (cm/s)
			ical.pca_CaMK = 1.1*ical.pca;				// Permiability of membrane to Ca (cm/s)
			ical.pcana_CaMK = 0.00125*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
			ical.pcak_CaMK = 3.574E-4*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
		} else if (var.celltype == 2){ // Mid.
			ical.pca = var.rategca*0.0001*2.5;      					// Permiability of membrane to Ca (cm/s)
			ical.pcana = 0.00125*ical.pca;   			// Permiability of membrane to Na (cm/s)
			ical.pcak = 3.574E-4*ical.pca;  			// Permiability of membrane to K (cm/s)
			ical.pca_CaMK = 1.1*ical.pca;				// Permiability of membrane to Ca (cm/s)
			ical.pcana_CaMK = 0.00125*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
			ical.pcak_CaMK = 3.574E-4*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
		} else if (var.celltype == 3){	// Epi.
			ical.pca = var.rategca*0.0001*1.2;      					// Permiability of membrane to Ca (cm/s)
			ical.pcana = 0.00125*ical.pca;   			// Permiability of membrane to Na (cm/s)
			ical.pcak = 3.574E-4*ical.pca;  			// Permiability of membrane to K (cm/s)
			ical.pca_CaMK = 1.1*ical.pca;				// Permiability of membrane to Ca (cm/s)
			ical.pcana_CaMK = 0.00125*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
			ical.pcak_CaMK = 3.574E-4*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
		}

	// Rapid delayed rectifier potassium current (Ikr)
		ikr.rategkr = sqrt(var.ko/5.4);

		if(var.celltype == 1){
			ikr.Gkr = 0.046;  //(mS/uF)
		} else if(var.celltype == 2){
			ikr.Gkr = 0.8*0.046;  //(mS/uF)
		} else if(var.celltype == 3){
			ikr.Gkr = 1.3*0.046;  //(mS/uF)
		}

	// Slow delayed rectifier potassium current (Ikr)
		if(var.celltype == 1){
			iks.Gks = 0.0034;  //(mS/uF)
		} else if(var.celltype ==2){
			iks.Gks = 0.0034;  //(mS/uF)
		} else if(var.celltype==3){
			iks.Gks = 1.4*0.0034;  //(mS/uF)
		}
	
	// Inward rectifier K current: Ik1
		if(var.celltype == 1){
			ik1.Gk1 = 0.1908; // (mS/uF)
		} else if(var.celltype == 2){
			ik1.Gk1 = 1.3*0.1908; // (mS/uF)
		} else if(var.celltype == 3){
			ik1.Gk1 = 1.2*0.1908; // (mS/uF)
		}
		ik1.rategk1 = sqrt(var.ko);

	// Sodium-Calcium Exchanger V-S
		var.kna1 = 15.0;		// (mM)
		var.kna2 =  5.0;		// (mM)
		var.kna3 = 88.12;		// (mM)
		var.kasym = 12.5;			// (mM)
		var.omega_na = 6.0E+4;		// (Hz)
		var.omega_ca = 6.0E+4;		// (Hz)
		var.omega_naca = 5.0E+3;	// (Hz)
		var.kca_on = 1.5E+6;		// (mM/ms)
		var.kca_off = 5.0E+3;		// (Hz)
		var.qna = 0.5224;   // 
		var.qca = 0.1670;   // 
		var.km_ca_act = 150.0E-6;		// (mM)

		if(var.celltype == 1){
			var.Gnaca = 0.0008;			// (uC/uF)
		} else if(var.celltype == 2){
			var.Gnaca = 1.4*0.0008;			// (uC/uF)
		} else if(var.celltype == 3){
			var.Gnaca = 1.1*0.0008;			// (uC/uF)
		}

	// Sodium-Potassium Pump
		inak.kp1 = 949.5;		// (Hz)
		inak.km1 = 182.4;		// (1/mM)
		inak.kp2 = 687.2;		// (Hz)
		inak.km2 = 39.4;			// (Hz)
		inak.kp3 = 1899.0;		// (Hz)
		inak.km3 = 79300.0;		// (Hz/mM/mM)
		inak.kp4 = 639.0;		// (Hz)
		inak.km4 = 40.0;			// (Hz)
		inak.ko_nai = 9.073;		// (mM)
		inak.ko_nao = 27.78;		// (mM)
		inak.delta = -0.1550;
		inak.kki = 0.5;			// (mM)
		inak.kko = 0.3582;		// (mM)
		inak.MgADP = 0.05;
		inak.MgATP = 9.8;
		inak.k_MgATP = 1.698E-7;	// (mM)
		inak.H = 1.0E-7;		// (mM) proton
		inak.SigP = 4.2;			// (mM)
		inak.HP = 1.698E-7;	// (mM)
		inak.nap = 224.0;		// (mM)
		inak.kp = 292.0;		// (mM)
		inak.b1 = inak.km1*inak.MgADP;
		inak.a2 = inak.kp2;
		inak.a4 = (inak.kp4*inak.MgATP/inak.k_MgATP)/(1.0 + inak.MgATP/inak.k_MgATP);

		if(var.celltype == 1){
			inak.G = 30.0;
		} else if(var.celltype == 2){
			inak.G = 0.7*30.0;
		} else if(var.celltype == 3){
			inak.G = 0.9*30.0;
		}

	// Sarcolemmal Ca Pump
		ipca.G = 0.0005;		// Max. Ca current through sarcolemmal Ca pump (mS/uF)
		ipca.km = 0.0005;		// Half-saturation concentration of sarcolemmal Ca pump (mM)

	// K Background Current 
		if(var.celltype == 1){
			ikb.G = 0.003;		// Max. conductance of K background (mS/uF)
		} else if(var.celltype==2){
			ikb.G = 0.003;		// Max. conductance of K background (mS/uF)
		} else if(var.celltype==3){
			ikb.G = 0.6*0.003;		// Max. conductance of K background (mS/uF)
		}

	// Ca Background Current 
		icab.pcab = 2.5E-8;		// (cm/s)
		icab.gacai = 1.0;
		icab.gacao = 0.341;

	// Na Background Current 
		inab.pnab = 3.75E-10;    // (cm/s)

	// Ca/CaM dependent protein kinese (CaMK)
		CaMK.a = 0.05;		// (1.0/ms)
		CaMK.b = 0.00068;	// (1.0/ms)
		CaMK.z = 0.05;
		CaM.Km = 0.0015;	// (mM)
		CaMK.Km = 0.15;
	
	// diffusion fluxes
		jdiff.tau_na = 2.0;	// (ms)
		jdiff.tau_k = 2.0;	// (ms)
		jdiff.tau_ca = 0.2;	// (ms)
	
	// SR calcium release flux, via RyR (Jrel)
		jrel.b_tau = 4.75;					// (ms)
		jrel.a = 0.5*jrel.b_tau; 			// (ms)
		jrel.b_tau_CaMK = 1.25*jrel.b_tau;	// (ms)
		jrel.a_CaMK = 0.5*jrel.b_tau_CaMK;	// (ms)
		if(var.celltype!=2){
			jrel.p = 1.0;
		} else {
			jrel.p = 1.7;
		}
	// calcium uptake via SERCA pump (Jup)
		jup.dKm_PLB = 0.00017;		// (mM)
		jup.dCaMK = 1.75;
		if(var.celltype == 1){
			jup.p = 1.0;
		} else if(var.celltype==2){
			jup.p = 1.0;
		} else if(var.celltype==3){
			jup.p = 1.3;
		}

	// Translocation of Ca Ions from NSR to JSR
		jtr.tau = 100.0;      // Time constant of Ca transfer from NSR to JSR (ms)

	// Myoplasmic Ca Ion Concentration Changes 
		if(var.celltype==1){
			var.cmdnbar = 0.050;   // Max. [Ca] buffered in CMDN (mM)
		} else if(var.celltype==2){
			var.cmdnbar = 0.050;   // Max. [Ca] buffered in CMDN (mM)
		} else if(var.celltype==3){
			var.cmdnbar = 1.3*0.050;   // Max. [Ca] buffered in CMDN (mM)
		}
		var.trpnbar = 0.070;   // Max. [Ca] buffered in TRPN (mM)
		var.bsrbar = 0.047;		// Max.[Ca] buffered in BSR (mM)
		var.bslbar = 1.124;		// Max.[Ca] buffered in BSL (mM)
		var.csqnbar = 10.0;     // Max. [Ca] buffered in CSQN (mM)
		var.kmcmdn = 0.00238;  // Equalibrium constant of buffering for CMDN (mM)
		var.kmtrpn = 0.0005;   // Equalibrium constant of buffering for TRPN (mM)
		var.kmbsr = 0.00087;  // Equalibrium constant of buffering for BSR (mM)
		var.kmbsl = 0.0087;   // Equalibrium constant of buffering for BSL (mM)
		var.kmcsqn = 0.8;     // Equalibrium constant of buffering for CSQN (mM)

}

