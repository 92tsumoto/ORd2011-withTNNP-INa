#include "syspara.h"

void eular(int n,double h,double x[],double t)
{
    //double k1[NN],k2[NN],k3[NN],k4[NN];
    double k1[NN];
    //double xtemp[NN];
    int i;  

    function(x,k1,t);

	//printf("%e ",t);
    for (i=0; i<n; i++){
		x[i]+=k1[i]*h;
	//	printf("%e ",x[i]);
	}
	//printf("\n");

}
