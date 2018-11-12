#include <string.h>
#include "syspara.h"

void eventloop(FILE *fp1,int *mode, int *P,double m[])
{
	int i;
		
	xhplot(KEYIN,0,0,WHITE);
	if (xddpret.key != -1){
		switch (xddpret.key){
			case 'e':
				if(*mode==1){
					xhplot(CLS, m[0],m[1],BLACK);
					xhplot(FRAME, m[0],m[1],WHITE);
				}
				if(*mode==2){
					xhplot(CLS, m[0],m[2],BLACK);
					xhplot(FRAME, m[0],m[2],WHITE);
				}
				if(!*mode){
					xhplot(CLS, m[0],m[1],BLACK);
					xhplot(FRAME, m[0],m[1],WHITE);
				}
				break;
			case 'A':
				var.Istim_base += var.dIstim;
				printf("Istim = %lf\n",var.Istim_base);
				fflush(stdout);
				break;
			case 'a':
				var.Istim_base -= var.dIstim;
				printf("Istim = %lf\n",var.Istim_base);
				fflush(stdout);
				break;
			case 'f':
				var.pflag = 1 - var.pflag;
				break;
			case 'g':
				switch(*P){
				case 2 : 
					*P = 8;
					fprintf(stderr,"fixed point is doted\n");
					printf("P=%d\n",*P);
					break;
				case 8 : 
					*P = 2;
					fprintf(stderr,"fixed point is not doted\n");
					printf("P=%d\n",*P);
					break;
				}
				break;
			case 'q':
				printf("quit\n");
				exit(0);
				break;
			case 's':
				printf("\nStatus\n");
				//printf("k2 = %lf --> v\n",var.k2);
				printf("%lf %lf %f %f %f\n", m[0],m[1],m[2],m[3],m[4]);
				printf("%lf %lf %f %f %f\n", m[5],m[6],m[7],m[8],m[9]);
				printf("%lf %lf %f %f %f\n", m[10],m[11],m[12],m[13],m[14]);
				printf("%lf %lf %f %f %f\n", m[15],m[16],m[17],m[18],m[19]);
				printf("%lf %lf %f %f %f\n", m[20],m[21],m[22],m[23],m[24]);
				printf("%lf %lf %f %f %f\n", m[25],m[26],m[27],m[28],m[29]);
				printf("%lf %lf %f %f %f\n", m[30],m[31],m[32],m[33],m[34]);
				printf("%lf %lf %e %e %f %f\n", m[35],m[36],m[37],m[38],m[39],m[40]);
				fflush(stdout);
				break;
			case '1':
				//fprintf(fp1,"%10.8lf\n",var.vd);
				//fprintf(fp1,"%10.8lf\n",var.k1);
				//fprintf(fp1,"%10.8lf\n",var.k2);
				for(i=0;i<NN;i++){fprintf(fp1,"%10.8lf\n",m[i]);}
				fflush(fp1);
				printf("stored in para.out\n");
				break; 
			case ' ':
				//mode = 1 - mode;
				printf("\n");
				*mode+=1;
				if(*mode==1)
					printf("### Now plane is Vm, m ###\n");
				if(*mode==2)
					printf("### Now plane is Vm, h ###\n");
				if(*mode>2){
					*mode=1; 
					printf("### Now plane is Vm, m ###\n");
					printf("\n");
					printf(" a:   Istim\n");
					printf(" e:   clear\n");
					printf(" f:   change\n");
					printf(" g:   small <-> big\n");
					printf(" s:   state\n");
					printf("space: Vm,m <-> Vm,h \n");
					printf("\n");
					break;
				}
		}
	}
}
