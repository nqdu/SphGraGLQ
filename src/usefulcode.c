#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*
    Convert geographic latitude to geocentric latitude
    hlat (input) = geographic latitude in radians (north positive)
    glat (output)= geocentric latitude in radians (north positive)
*/
double geogra2geocenlat(double hlat)
{
    double halfpi = 1.570796327;
    double x;

    if ((halfpi - fabs(hlat)) >= 0.05) {    //===far away from the pole==
        //x = atan(0.993277*sin(hlat)/cos(hlat));
        x = atan(0.993277*tan(hlat));
    }
    else {
        /* Special formula near pole */
        if (hlat > 0)
        x  = hlat/0.993277 - 0.010632;
        else
        x  = hlat/0.993277 + 0.010632;
    }

    return x;
}

//------convert geocentric latitude back to geographic latitude-------------
//       glatinv (output) = geographic latitude in radians (north positive)
//       hlat (input)= geocentric latitude in radians (north positive)
//---------------------------------------------------------------------
double  geocen2geogralat(double hlat)
{
    double halfpi = 1.570796327;
    double glatinv;
    if( (halfpi - fabs(hlat)) >= 0.05){
        //glatinv = atan(sin(hlat)/cos(hlat)/0.993277);
        glatinv = atan(tan(hlat)/0.993277);
    }
    else{   // ------special formula near pole
        if( hlat>0 )
            glatinv = (hlat+0.010632)*0.993277;
        else
            glatinv = (hlat-0.010632)*0.993277;
    }
    return glatinv ;
}

/*
    compute straight line distance between points x_e and x_s points
    the lon,colat,r in rad and km are used as inputs
*/
double compute_length(double colatErad, double lonErad, double RE, 
                    double colatSrad, double lonSrad, double RS )
{
    double R = 6371.0;
    double PIhalf = 1.5707963;
    double Lates, Lones, Res, rep, d2;
    double latErad, latSrad;
    double length;

    //=====Modified 2010/02/15
    latErad = PIhalf - colatErad;
    latSrad = PIhalf - colatSrad;

    Lates   = colatErad - colatSrad;
    Lones   = lonErad - lonSrad;
    Res     = RE - RS;

    rep = Lones * sin( (latErad+latSrad)/2.0 );

    //d2  = Lates*Lates+rep*rep;
    //(*length) = sqrt( Res*Res+d2*RE*RS );
    length = sqrt( Res*Res + (Lates*Lates+rep*rep)*RE*RS );

    return length;
}

/*
	calculate the minmum distance to one grid surface
*/
void cal_minDis(double *mindis, double colatrado, double lonrado, double ro, 
            double colatrads1, double lonrads1, double colatrads2, double lonrads2, 
            double rs )
{
	double length[4];

	/*
		if observation point is outside the same lon/lat
	*/
	if( (colatrado<colatrads1||colatrado>colatrads2) && (lonrado<lonrads1||lonrado>lonrads2) ){
		length[0]=compute_length(colatrado, lonrado, ro, colatrads1, lonrads1, rs ); 
		length[1]=compute_length(colatrado, lonrado, ro, colatrads1, lonrads2, rs );
		length[2]=compute_length(colatrado, lonrado, ro, colatrads2, lonrads1, rs );
		length[3]=compute_length(colatrado, lonrado, ro, colatrads2, lonrads2, rs );

		(*mindis) = length[0];
        if( (*mindis) > length[1] ) (*mindis) = length[1];
        if( (*mindis) > length[2] ) (*mindis) = length[2];
        if( (*mindis) > length[3] ) (*mindis) = length[3];
	}
	/*
        if the only longitude of observation point is inside the same lon
    */
	else if((lonrado>=lonrads1&&lonrado<=lonrads2) && (colatrado<colatrads1||colatrado>colatrads2)  ){
		if( fabs(colatrado-colatrads1)<fabs(colatrado-colatrads2) ){
			(*mindis)=compute_length(colatrado, lonrado, ro, colatrads1, lonrado, rs );
		}
		else{
			(*mindis)=compute_length(colatrado, lonrado, ro, colatrads2, lonrado, rs );
		}

	}
	/*
                if the only latitude of observation point is inside the same lat
        */
	else if((colatrado>=colatrads1&&colatrado<=colatrads2) && (lonrado<lonrads1||lonrado>lonrads2) ){
		if( fabs(lonrado-lonrads1)<fabs(lonrado-lonrads2) ){
			(*mindis)=compute_length(colatrado, lonrado, ro, colatrado, lonrads1, rs );
		}
		else{
			(*mindis)=compute_length(colatrado, lonrado, ro, colatrado, lonrads2, rs );
		}
	}
	else{
        (*mindis) = fabs(ro-rs);
	}

}

double *dvec(int n){
    double *a=(double*)calloc(sizeof(double),n);
    return a;
}

float *fvec(int n){
    float *a=(float *)calloc(sizeof(float),n);
    return a;
}

int *ivec(int n){
    int *a=(int *)calloc(sizeof(int),n);
    return a;
}


float **fmat(int row,int col){
    int i;
    float **a=(float **)calloc(sizeof(float*),row);

    for(i=0;i<row;i++) a[i]=fvec(col);
    return a;
}

int **imat(int row,int col){
    int i;
    int **a=(int **)calloc(sizeof(int*),row);

    for(i=0;i<row;i++) a[i]=ivec(col);

    return a;
}

double **dmat(int row,int col){
    int i;
    double **a=(double **)calloc(sizeof(double*),row);

    for(i=0;i<row;i++) a[i]=dvec(col);
    return a;
}

double ***d3tensor(int row,int col,int deep){
    int i;
    double ***a=(double ***)calloc(sizeof(double**),row);

    for(i=0;i<row;i++) a[i]=dmat(col,deep);
    return a;
}

void free_dvec(double *a){
    free(a);
    a=NULL;
}

void free_fvec(float *a){
    free(a);
    a=NULL;
}

void free_ivec(int *a){
    free(a);
    a=NULL;
}

void free_dmat(double **P,int row_p){
    int i;
    for(i=0;i<row_p;i++) free_dvec(P[i]);
    free(P);
    P=NULL;
}

void free_fmat(float **a,int row){
    int i;
    for(i=0;i<row;i++) free_fvec(a[i]);
    free(a);
    a=NULL;
}

void free_imat(int **a,int row){
    int i;
    for(i=0;i<row;i++) free_ivec(a[i]);
    free(a);
    a=NULL;  
}

void free_d3tensor(double ***P,int row,int col){
    int i;
    for(i=0;i<row;i++) free_dmat(P[i],col);
    free(P);
    P=NULL;
}
