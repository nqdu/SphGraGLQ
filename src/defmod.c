#include <stdlib.h>
#include<stdio.h>
#include<string.h>
#include"defmod.h"
#include"usefulcode.h"

void free_velomodel(MOD3DSphGra *mod3d)
{
    free_fvec(mod3d->lon);
    free_fvec(mod3d->lat);
    free_fvec(mod3d->dep);
    free_fvec(mod3d->density0);
}

void free_OBSSphGraRandom(OBSSphGraRandom *obs)
{
    free_fvec(obs->lon);
    free_fvec(obs->lat);
    free_fvec(obs->Gr);
}

// input and output function
void chancoorobssphonlygra( OBSSphOnlyGra *ObsSphOnlyGra,int israd )
{
    int     i, nx, ny, nz;
	double  lat0Geocen, dlatrad;

    nx = ObsSphOnlyGra->nx;
    ny = ObsSphOnlyGra->ny;

    if( israd==0 && ObsSphOnlyGra->israd==1 ){          //change coordinates to lon/lat/dep;
        for(i=0; i<nx; i++ ){
                ObsSphOnlyGra->lon[i] = ObsSphOnlyGra->lon[i]/degrad2;
        }


        for(i=0; i<ny; i++ )
            ObsSphOnlyGra->lat[i] =  geocen2geogralat((double)(hpi2 -
                                        ObsSphOnlyGra->lat[i]) )/degrad2;

        ObsSphOnlyGra->x0   = ObsSphOnlyGra->x0/degrad2;
        //ObsSphOnlyGra->y0   = (hpi2 - ObsSphOnlyGra->y0)/degrad2;
        ObsSphOnlyGra->y0   = geocen2geogralat((double)(hpi2 - ObsSphOnlyGra->y0))/degrad2;

        ObsSphOnlyGra->z0   = rearth2 - ObsSphOnlyGra->z0;

        //ObsSphOnlyGra->dlon = ObsSphOnlyGra->dlon/degrad2;
        //ObsSphOnlyGra->dlat = ObsSphOnlyGra->dlat/degrad2;

        ObsSphOnlyGra->israd = 0;
    }
	else if( israd==1 && ObsSphOnlyGra->israd==0) {     //change coordinates to lonrad/colatrad/r;
        for(i=0; i<nx; i++ )
            ObsSphOnlyGra->lon[i] = ObsSphOnlyGra->lon[i]*degrad2;

        for(i=0; i<ny; i++ ){
            ObsSphOnlyGra->lat[i] = hpi2 - geogra2geocenlat( (double)(ObsSphOnlyGra->lat[i]*degrad2) );
        }

        ObsSphOnlyGra->x0   = ObsSphOnlyGra->x0*degrad2;
        //ObsSphOnlyGra->y0   = hpi2 - ObsSphOnlyGra->y0*degrad2;
        ObsSphOnlyGra->y0   = hpi2 - geogra2geocenlat((double)(ObsSphOnlyGra->y0*degrad2) );
        ObsSphOnlyGra->z0   = rearth2 - ObsSphOnlyGra->z0;

        ObsSphOnlyGra->israd = 1;
    }
    else{
        printf( "israd=%d and ObsSphOnlyGra->israd=%d\n", israd, ObsSphOnlyGra->israd );
        printf( "Needn't change coordinates.\n");
    }
    
}

void chancoorsphgramod(MOD3DSphGra *mod3dsphgra, int israd )

{
	int	i, nx, ny, nz;

	nx = mod3dsphgra->nx;
    ny = mod3dsphgra->ny;
    nz = mod3dsphgra->nz;    

	if( israd==0 && mod3dsphgra->israd==1 ){//change coordinates to lon/lat/dep;
		for(i=0; i<nx; i++ )
            mod3dsphgra->lon[i] = mod3dsphgra->lon[i]/degrad2;

        for(i=0; i<ny; i++ )
			mod3dsphgra->lat[i] = geocen2geogralat((double)(hpi2 - mod3dsphgra->lat[i]) )/degrad2;

        for(i=0; i<nz; i++)
            mod3dsphgra->dep[i] = rearth2 - mod3dsphgra->dep[i];

		mod3dsphgra->x0	  = mod3dsphgra->x0/degrad2;
		//mod3dsphgra->y0   = (hpi2 - mod3dsphgra->y0)/degrad2;
		mod3dsphgra->y0   = geocen2geogralat((double)(hpi2 - mod3dsphgra->y0) )/degrad2;

		mod3dsphgra->z0   = rearth2 - mod3dsphgra->z0;

		mod3dsphgra->israd = 0;
	}
	else if( israd==1 && mod3dsphgra->israd==0) {//change coordinates to lonrad/colatrad/r;
		for(i=0; i<nx; i++ )
            mod3dsphgra->lon[i] = mod3dsphgra->lon[i]*degrad2;

        for(i=0; i<ny; i++ )
			mod3dsphgra->lat[i] =   hpi2 - geogra2geocenlat( (double)(mod3dsphgra->lat[i]*degrad2) );


        for(i=0; i<nz; i++)
            mod3dsphgra->dep[i] = rearth2 - mod3dsphgra->dep[i];

		mod3dsphgra->x0   = mod3dsphgra->x0*degrad2;
		mod3dsphgra->y0   = hpi2 - geogra2geocenlat( (double)(mod3dsphgra->y0*degrad2) );
        mod3dsphgra->z0   = rearth2 - mod3dsphgra->z0;

		mod3dsphgra->israd = 1;
	}
	else{
		printf( "israd=%d and mod3dsphgra->israd=%d\n", israd, mod3dsphgra->israd );
		printf( "Needn't change coordinates.\n");
	}

}

/*
        if ISRAD=0, change coordinates to lon/lat/dep;
        if ISRAD=1, change coordinates to lonrad/colatrad/r;
*/
void chancoorOBSsphGraRandom(OBSSphGraRandom *ObsSphGra, int israd )
{
    int     i, np, nz;
    np = ObsSphGra->np;

    if( israd==0 && ObsSphGra->israd==1 ){          //change coordinates to lon/lat/dep;
        for(i=0; i<np; i++ ){
            ObsSphGra->lon[i] = ObsSphGra->lon[i]/degrad2;
            //ObsSphGra->lat[i] = (hpi2 - ObsSphGra->lat[i])/degrad2;
            ObsSphGra->lat[i] =  geocen2geogralat((double)(hpi2 - ObsSphGra->lat[i]) )/degrad2;
        }
        ObsSphGra->z0   = rearth2 - ObsSphGra->z0;
        ObsSphGra->israd = 0;
    }
    else if(israd==1 && ObsSphGra->israd==0) {     //change coordinates to lonrad/colatrad/r;
        for(i=0; i<np; i++ ){
            ObsSphGra->lon[i] = ObsSphGra->lon[i]*degrad2;
            //ObsSphGra->lat[i] = hpi2 - ObsSphGra->lat[i]*degrad2;
            ObsSphGra->lat[i] = hpi2 - geogra2geocenlat( (double)(ObsSphGra->lat[i]*degrad2));
        }
        ObsSphGra->z0   = rearth2 - ObsSphGra->z0;
        ObsSphGra->israd = 1;
    }
	else{
            fprintf(stderr, "ISRAD=%d and ObsSphGra->ISRAD=%d\n", israd, ObsSphGra->israd );
            fprintf(stderr, "Needn't change coordinates.\n");
    }
}

/*
    Read observed gravity Gr data, the data file should be
    in 3 columns like (lon,lat,Gr),and the number of obeservation
    is in the front of this file
*/
void readOBSSphGraRandom(OBSSphGraRandom *obs,char *filename)
{
    FILE *fp;
    int i,np;
    fp=fopen(filename,"r");

    fscanf(fp,"%d",&np); // read number of observations
    obs->np=np;
    obs->z0=0.0;

    // allocate space
    obs->lat = fvec(np);
    obs->lon= fvec(np);
    obs->Gr = fvec(np);

    // read obervations
    for(i=0;i<np;i++)
        fscanf(fp,"%f%f%f",obs->lon+i,obs->lat+i,obs->Gr+i);
    
    fclose(fp);
    obs->israd=0;
}

// read model parameters and allocate space
// note that the initial model is larger than inverted one
int read_MOD3DSphGra(MOD3DSphGra *mod3d,char *paramfile,char *modinfile)
{   
    FILE *fp,*fptrue;
    int i,j,k;
    int nx,ny,nz,n,kmax;
    float ulx,uly,dx,dy;
    float v,vtrue,maxdis;
    char line[256];

    fp=fopen(paramfile,"r");
    for(i=0;i<4;i++)
        fgets(line,256*sizeof(char),fp);
    fgets(line,256*sizeof(char),fp);

    sscanf(line,"%d%d%d",&ny,&nx,&nz);
    mod3d->nx=nx-2;
    mod3d->ny=ny-2;
    mod3d->nz=nz-1;
    n=(nx-2)*(ny-2)*(nz-1);
    mod3d->lon = fvec(nx-2);
    mod3d->lat = fvec(ny-2);
    mod3d->dep = fvec(nz-1);

    fgets(line,256*sizeof(char),fp);
    sscanf(line,"%f%f",&uly,&ulx);
    fgets(line,256*sizeof(char),fp);
    sscanf(line,"%f%f",&dy,&dx);
    ulx=ulx+dx;
    uly=uly-dy;
    mod3d->x0=ulx;
    mod3d->y0=uly;

    for(i=0;i<nx-2;i++)
        mod3d->lon[i]=ulx+i*dx;
    for(i=0;i<ny-2;i++)
        mod3d->lat[i]=uly-i*dy;

    //read synflag
    int synflag;
    for(i=0;i<7;i++)
        fgets(line,256*sizeof(char),fp);
    sscanf(line,"%d",&kmax); //kmaxRc
    if(kmax>0){
        fgets(line,256*sizeof(char),fp);
        fgets(line,256*sizeof(char),fp);
    }
    else{
        fgets(line,256*sizeof(char),fp);
    }

    sscanf(line,"%d",&kmax); //kmaxRg
     if(kmax>0){
        fgets(line,256*sizeof(char),fp);
        fgets(line,256*sizeof(char),fp);
    }
    else{
        fgets(line,256*sizeof(char),fp);
    } 

    sscanf(line,"%d",&kmax); //kmaxLc
     if(kmax>0){
        fgets(line,256*sizeof(char),fp);
        fgets(line,256*sizeof(char),fp);
    }
    else{
        fgets(line,256*sizeof(char),fp);
    }     

    sscanf(line,"%d",&kmax); //kmaxLg
     if(kmax>0){
        fgets(line,256*sizeof(char),fp);
        fgets(line,256*sizeof(char),fp);
    }
    else{
        fgets(line,256*sizeof(char),fp);
    }   

    sscanf(line,"%d",&synflag);
    //read maxdis for calculation of GLQ points    
    for(i=0;i<3;i++)
        fgets(line,256*sizeof(char),fp);
    sscanf(line,"%f",&maxdis);
    fclose(fp);

    /*---------------------------------------------------*/
    // read depth of velomodel
    fp=fopen(modinfile,"r");
    for(i=0;i<nz-1;i++)
        fscanf(fp,"%f",mod3d->dep+i);
    mod3d->density0 = fvec(n);

    fscanf(fp,"%f",&v); // read dummy depth
    if(synflag==1){
        mod3d->density = fvec(n);
        strcat(modinfile,".true");
        fptrue=fopen(modinfile,"r");
    }
    
    // velocity model has shape (dep,lon,lat)
    for(k = 0;k < nz;k++){
        for(i = 0;i < nx;i++){
            for(j=0;j<ny;j++){
                n = k*(nx-2)*(ny-2)+(i-1)*(ny-2)+(j-1);
                fscanf(fp,"%f",&v);
                if(synflag==1)
                    fscanf(fptrue,"%f",&vtrue);
                if(i == 0||j == 0||k == nz-1||i == nx-1||j == ny-1)
                    continue;
                mod3d->density0[n]=v*1.732*0.3601+0.541;
                if(synflag==1)
                    mod3d->density[n]=vtrue*1.732*0.3601+0.541;
            }
        }
    }
    if(synflag == 1 )
        fclose(fptrue);
    mod3d->israd = 0;
    fclose(fp);

    return synflag;

}