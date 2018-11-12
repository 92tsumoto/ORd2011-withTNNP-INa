/* produced by Tsumoto. K 2008.10.27 */

#include <string.h>
#include <stdlib.h>
#include "syspara.h"

FILE *fopen(), *fpin, *fp0, *fp1, *fp2, *fp3;
int mode = 1;
int P = 2;
int beats = 30;

typedef double Number;

main(argc,argv)
int argc;
char **argv;
{
	int i,w;
	int ii=0;
	double x[NN];
	double t = 0.0;
	double time=0.0;
	double h;
	double v_old,dvdt,dvdt_new;
	double t_stok;
	char *tmpname;
	char cmd[BUFSIZ];
	double tend;

/* Action Potential Duration and Max. Info */
	double *vmax ; // Max. Voltage (mV)
	double *dvdtmax ; // Max. dv/dt (mV/ms)
	double *apd; // Action Potential Duration
	double *toneapd; // Time of dv/dt Max.
	double *ttwoapd; // Time of 90% Repolarization
	double *rmbp; // Resting Membrane Potential
	double *nair; // Intracellular Na At Rest
	double *cair; // Intracellular Ca At Rest
	double *kir ; // Intracellular K At Rest
	double caimax [beats] ; // Peak Intracellular Ca

	vmax=(Number *)calloc(beats,sizeof(Number));
	dvdtmax=(Number *)calloc(beats,sizeof(Number));
	apd=(Number *)calloc(beats,sizeof(Number));
	toneapd=(Number *)calloc(beats,sizeof(Number));
	ttwoapd=(Number *)calloc(beats,sizeof(Number));
	rmbp=(Number *)calloc(beats,sizeof(Number));
	nair=(Number *)calloc(beats,sizeof(Number));
	cair=(Number *)calloc(beats,sizeof(Number));
	kir=(Number *)calloc(beats,sizeof(Number));
	if(vmax==NULL || dvdtmax==NULL || apd==NULL || toneapd==NULL || ttwoapd==NULL 
		|| rmbp==NULL || nair==NULL || cair==NULL || kir==NULL
		) exit(1);

//int i; // Stimulation Counter

	tmpname = "temp";

	sprintf(cmd, "/usr/bin/cpp -P %s > %s", argv[1],tmpname);
	if(system(cmd) == -1){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if((fpin=fopen(tmpname,"r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if ((fp1 = fopen("para.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp2 = fopen("data.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp3 = fopen("ndata.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}

// parameter inputs
	input_para(fpin);

	if (var.write){
		if ((fp0 = fopen(argv[2],"w"))==NULL){
			fprintf(stderr, "%s cannot open.\n",argv[2]);
			exit(-1);
		}
	}

	xhplot(WINDOW, 700.0, 700.0, WHITE);
	xhplot(DIRECT, 0.0, 0.0, WHITE);

	for (ii = 0; ii < var.datas; ii++){
		long j;
		time = 0.0;
		tend = var.tend[ii];
		for (i = 0; i < NN; i++){ 
			x[i] = var.x0[ii][i];
		}
		if (var.half){
			h = 1.0 / var.m;
		}
		else {
			h = 1.0 / var.m;
		}
		h *= var.tsign[ii];
		xddp.line_wid = var.line_wid[ii];
		xhplot(LINEATT,0,0,WHITE);

		
		// initial values input.
		val_consts(fp1);
		printf("exit consts\n");
	
	// initial values input.
		initial_mem();
		printf("exit memory initialization\n");

		printf("Istim=%lf\n",var.Istim_base);

	// Tablize exp functions.	
	printf("start tablization\n");
	make_ExpTable();
	printf("finished tablization\n");

	// Initialization time
		time -= h;
		var.dt = h;
		var.beat = 0;

		ii = 0;
		
		while (1){
			eventloop(fp1,&mode,&P,x);

			for (j = 0; j< (var.m * var.l ); j++){
				t = h*j;
				time += h;

				if ( time-(var.BCL*var.beat+10.0) >= 0.0 && time-(var.BCL*var.beat+10.0) < h ){
					apd[var.beat] =0;
					toneapd[var.beat] =0;
					ttwoapd[var.beat] =0;
					rmbp[var.beat] =x[0];
					nair[var.beat] = x[32];
					kir[var.beat] = x[34];
					cair[var.beat] = x[36];
					caimax[var.beat] = x[36];
					vmax[var.beat] = -90.0;
					dvdtmax[var.beat] = 0.0;

					printf("%d %lf %lf %lf %lf\n",var.beat,apd[var.beat],toneapd[var.beat],ttwoapd[var.beat],rmbp[var.beat]);
					printf("%lf %lf %f %f %e\n", x[0],x[1],x[2],x[3],x[4]);
					printf("%lf %lf %f %f %e\n", x[5],x[6],x[7],x[8],x[9]);
					printf("%lf %e %lf %lf %lf\n", x[10],x[11],x[12],x[13],x[14]);
					printf("%lf %e %lf %lf %lf\n", x[15],x[16],x[17],x[18],x[19]);
					printf("%lf %e %lf %lf %lf\n", x[20],x[21],x[22],x[23],x[24]);
					printf("%lf %e %lf %lf %lf\n", x[25],x[26],x[27],x[28],x[29]);
					printf("%lf %e %e %lf %lf\n", x[30],x[31],x[32],x[33],x[34]);
					printf("%lf %lf %e %e %e\n", x[35],x[36],x[37],x[38],x[39]);
					printf("time=%lf,Istim=%lf\n",time,var.Istim);
					//printf("dvdtmax[%d]=%lf\n",var.beat,dvdtmax[var.beat]);
					printf("dvdtmax[%d]=%lf\n",var.beat,dvdt);
					printf("ENa=%lf, EK=%lf, EKs=%lf\n", var.Ena,var.Ek,var.Eks);
				}

				//if(var.beat < 30){
				if (time-(var.BCL*var.beat+10.0) >= 0.0 && time-(var.BCL*var.beat+10.0) < 0.5){
					var.Istim = var.Istim_base;
				} else {
					var.Istim = 0;
				}
				//}	

				if (fabs(time) > tend &&  tend != 0.0) break;

				v_old = x[0];

				eular(NN,h,x,t);
				
			/*	printf("%e ",time);
				for (i=0; i<NN; i++){
				printf("%e ",x[i]);
				}printf("%e\n",var.Istim);
			*/
				dvdt_new = (x[0]-v_old)/h;
				//printf("dvdt_new=%lf\n",dvdt_new);

				if(var.beat>=0){
					if (x[0] > vmax[var.beat] )
						vmax[var.beat] = x[0];
					if (x[36] > caimax[var.beat] )
						caimax[var.beat] = x[36];
					if (dvdt_new > dvdtmax[var.beat] ){
						dvdtmax[var.beat] = dvdt_new;
						toneapd[var.beat] = time;
					}
					if (dvdt_new < 0 && x[0] >= (vmax[var.beat] -0.9*(vmax[var.beat]-rmbp[var.beat]) ) )
						ttwoapd[var.beat] = time;
				}

				if (var.pflag) orbit(&mode,x,dvdt_new);

				if (time>= (beats-5)*var.BCL && time < beats*var.BCL){
				//if (time>= 0.0){
					fprintf(fp2,"%lf %lf\n",time,x[0]);
					//fprintf(fp2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %e %e %e %e %e %e\n",
					//			time,x[0],ina.fast,ina.late,ito.ik,ikr.ik,iks.ik,ik1.ik,var.inaca,var.inaca_ss,ical.ica,inak.inak,inab.na,icab.ca,x[33],x[34],x[37],x[38],x[39],x[40]);
					//fprintf(fp2,"%lf %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
					//			time,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],
					//			x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],
					//			x[21],x[22],x[23],x[24],x[25],x[26],x[27],x[28],x[29],x[30],
					//			x[31],x[32],x[33],x[34],x[35],x[36],x[37],x[38],x[39],x[40],var.Istim);
				}
				
				dvdt = dvdt_new;

			}
			ii++;

			if(ii==var.BCL){
				var.beat++;
				ii = 0;
				if(var.beat > beats){
					for(w=0;w<=beats;w++){
					}
					for(w=0;w<=beats;w++){
						apd[w] = ttwoapd [w] -toneapd [w] ;
						fprintf(fp3,"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%e\t%e\t%g\n",w,
							vmax[w],dvdtmax[w],apd[w],toneapd[w],ttwoapd[w],nair[w],kir[w],cair[w],caimax[w],rmbp[w]);
						printf("%d %lf %lf %lf %lf %lf %lf %lf %e %e %lf\n",w,
							vmax[w],dvdtmax[w],apd[w],toneapd[w],ttwoapd[w],nair[w],kir[w],cair[w],caimax[w],rmbp[w]);
					}
					exit(0);
				}
			}

			draw_p(&mode,P,x,dvdt);
			mouse(&mode,x,dvdt);
			if (fabs(time) > tend &&  tend != 0.0) break;

		}
		fclose(fp1);
		fclose(fp2);
		fclose(fp3);
		free(vmax);free(dvdtmax);free(apd);free(toneapd);free(ttwoapd);
		free(rmbp);free(nair);free(cair);free(kir);
		closed_mem();
	}
}

