#include"defmod.h"

void OBSGravityGrRandomAdaptive(MOD3DSphGra *mod3dsphgra,OBSSphGraRandom *ObsSphGra );
int gravmat(MOD3DSphGra *mod3dsphgra,OBSSphGraRandom *ObsSphGra,float *rw,int *iw);
int gravmat_parallel(MOD3DSphGra *mod3dsphgra,OBSSphGraRandom *ObsSphGra,float *rw,int *iw);