#include<stdio.h>
#include<stdlib.h>
#include"forward.h"
#include"usefulcode.h"
#include<math.h>
/*
	compute y=y+A*x,
	m,n   : shape of matrix A
	leniw : length of iw
	lenrw : length of rw
	iw	  : iw[0] is the number of nonzeros in A,
			iw[1:iw[0]] is corresponding row,
			iw[iw[0]+1,:2*iw[0]+1] is corresponding col
	rw    : non-zero elements of A 
*/
void aprod(int m, int n, float *x, float *y, 
			int leniw, int lenrw, int *iw, float *rw)
{
	int k;
	int nar=iw[0];

	for(k=0;k<nar;k++){
		y[iw[k+1]]+=rw[k]*x[iw[nar+1+k]];
	}
}

int main(){

	MOD3DSphGra mod3d;
	OBSSphGraRandom obsgrav;
	int i,j,k;
	int synflag;                                                                                                                                                                                                        
	int nmd,ndat;
	char modinfile[256], paramfile[256],obsinfile[256];
	FILE *fp;

	sprintf(modinfile,"MOD");
	sprintf(paramfile,"DSurfTomo.in");
	sprintf(obsinfile,"obsgrav.dat");

	// read model parameters and allocate space
	// note that the initial model is larger than inverted one
	synflag = read_MOD3DSphGra(&mod3d,paramfile,modinfile);
	
	// read observation data
	readOBSSphGraRandom(&obsgrav,obsinfile);
	printf("finish reading\n");

	// change degrees to radians
	chancoorOBSsphGraRandom(&obsgrav,1);
	chancoorsphgramod(&mod3d,1);

	// generate gravity matrix
	if(synflag == 1)
		printf("begin to compute gravity matrix and also synthesize gravity data\n");
	else
		printf("begin to compute gravity matrix\n");
	
	//OBSGravityGrRandomAdaptive(&mod3d,&obsgrav);
	nmd = mod3d.nx*mod3d.ny*mod3d.nz;
	ndat = obsgrav.np;
	float *rw;
	int *iw;
	rw = (float *)calloc(sizeof(float),nmd*ndat);
	rw = fvec(nmd*ndat);
	iw= ivec(2*nmd*ndat+1);
	float *y = fvec(ndat);
	float *x = fvec(nmd);
	
	int nar = gravmat_parallel(&mod3d,&obsgrav,rw,iw);
	if(synflag == 1){
		for(i = 0;i < nmd;i++)
		   x[i] = mod3d.density[i] - mod3d.density0[i];
		aprod(ndat,nmd,x,y,2*nar+1,nar,iw,rw);
	}
	else{
		for(i = 0;i < ndat;i++)
			y[i] = obsgrav.Gr[i];
	}
	// output
	chancoorOBSsphGraRandom(&obsgrav,0);
	fp=fopen("gravity.dat","w");
	fprintf(fp,"%d\n",ndat);
	for(i=0;i<obsgrav.np;i++)
		fprintf(fp,"%f %f %g\n",obsgrav.lon[i],obsgrav.lat[i],y[i]);
	fclose(fp);
	
	fp=fopen("gravmat.dat","w");  
	fprintf(fp,"%d\n",nar);
	for(i=0;i<nar;i++)
		fprintf(fp,"%d %d %g\n",iw[i+1],iw[i+nar+1],rw[i]);
	fclose(fp);

	// deallocate space
	free_OBSSphGraRandom(&obsgrav);
	free_velomodel(&mod3d);
	if(synflag==1){
		free_fvec(mod3d.density);
	}

	free_fvec(rw);
	free_ivec(iw);
	free_fvec(y);
	free_fvec(x);

	return 0;
}
