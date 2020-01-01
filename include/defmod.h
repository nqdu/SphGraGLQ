#include<stdlib.h>
#define  degrad2 0.017453293
#define  hpi2    1.570796327
#define  rearth2 6371.0
#define G 6.67E-6
#define pi 3.141592654

//=== define model and data structure ===

/*	This model is used to compute gravity effects and compute traveltime*/
typedef struct{
    int     nx, ny, nz;            		//=== number of inversion grids in x, y, z, directions ===
    float   x0, y0, z0;
    float   *lon, *lat, *dep;
	float   *density0, *density;

	int	israd;				//===Flag for lon/lat/dep (israd=0) or lonrad/colatrad/r (israd=1)====
}MOD3DSphGra;

//=== Uniform observation system, all variables are in Observation-Centred Coordinate System
typedef struct{
    int     nx, ny;                         //=== number of observation grids in x, y, directions ===
    float   x0, y0;
    float   z0;                             //=== Depth, Downward is positive; NOTICE
    float   dlon, dlat;                     //=== intervals in x, y direction
    float   *lon, *lat;

    float  *gravity;         		//=== gravity Vector components in radialdirection

    int     israd;   			//=== Flag for lon/lat (israd=0) or lonrad/colatrad/r (israd=1)====
}OBSSphOnlyGra;

//=== Random observation system, all variables are in Observation-Centred Coordinate System
typedef struct{
    int     np;        		//=== number of observation grids in x, y, directions ===
	float	z0;			//=== observation elevation
    float   *lon, *lat;		//=== observation lon/lat

	float	*Gr;					//=== Observered gravtiy in radial direction.

    int     israd;                                  //=== Flag for lon/lat (ISRAD=0) or lonrad/colatrad/r (ISRAD=1)====
}OBSSphGraRandom;

//===========for lsqr algorithm====
typedef struct
{
    long int        npt; 		// number of none-zero elements. /
    long int        *idx;		// array to store the index of none-zero value
    float          *value; 	// array to store the none-zero value & idx

}Srow;

typedef struct
{
    long int        nrow, ncol;     // row number & colume number
    Srow            *smx;

}Smat;

//deallocate memory function
void free_velomodel(MOD3DSphGra *mod3d);
void free_obs(OBSSphOnlyGra *obs);
void free_Smat(Smat *MA);
void free_OBSSphGraRandom(OBSSphGraRandom *obs);

// read observation gravity function
void readOBSSphOnlyGraBinary( OBSSphOnlyGra *ObsSphOnlyGra, char *filename );
void readOBSSphGraRandom(OBSSphGraRandom *obs,char *filename);
int read_MOD3DSphGra(MOD3DSphGra *mod3d,char *paramfile,char *modinfile);

// change degree to rad
void chancoorobssphonlygra( OBSSphOnlyGra *ObsSphOnlyGra,int israd );
void chancoorsphgramod(MOD3DSphGra *mod3dsphgra, int israd );
void chancoorOBSsphGraRandom(OBSSphGraRandom *ObsSphGra, int israd );
