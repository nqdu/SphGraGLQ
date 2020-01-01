#include<stdio.h>
#include<math.h>
#include"defmod.h"
#include"gauleg.h"
#include"usefulcode.h"
/*
	Change NMAX and NSCALE bigger for more accurate results,
	but bigger NMAX and NSCALE means much more computation.
*/
//===minimum and maximum nodes for GLQ
//===the ratio of distance/GLQ grid spacing
//===the minmum scale for GLQ, in km
#define	NMIN      2
#define NMAX      256
#define NSCALE    5
#define MINDISKM  0.5

/*
    Given position in lonrad/colatrad/r, get the density or the anomally density
*/
void getdensityBlock(MOD3DSphGra *mod3dsphgra, int px, int py, int pz, double *density )
{
    int     nx, ny, nz;
    int     i, j, k;
    long    N, nxy;
                                                                                                                                                                                
    nx = mod3dsphgra->nx;
    ny = mod3dsphgra->ny;
    nz = mod3dsphgra->nz;
    nxy = nx*ny;
                                                                                                                                                                                   
    if( px<0 || px>=nx || py<0 || py>=ny || pz<0 || pz>=nz ){
            fprintf(stderr,"point out of model. px, py, pz=%d %d %d\n", px, py, pz );
            fprintf(stderr,"model:nx, ny, nz = %d %d %d\n", nx, ny, nz );
            exit(0);
    }
                                                                                                                                                                                
    N = pz * nxy + px * ny + py;
    (*density) = (double)(mod3dsphgra->density[N]-mod3dsphgra->density0[N]);
}

/*
	for Radial Gravity, from GravityEffects
*/
void GravityR(double ro,double lonrado,double colatrado,
            double rs,double lonrads,double colatrads, 
            double density,double *Gr)
{       
        double  K, costmp, Ros, Ros3;
        
        //====change to International Unit kg/m3
        density = density*1000.0;
        ro      = ro*1000.0;
        rs      = rs*1000.0;
        
        K = G*density*rs*rs*sin(colatrads);
        costmp = cos(colatrado)*cos(colatrads) + sin(colatrado)*sin(colatrads)*cos(lonrado-lonrads );
        Ros    = sqrt(ro*ro + rs*rs - 2.0*ro*rs*costmp );
        
        Ros3 = Ros*Ros*Ros;
        
         //NOTE: According the Observation coordinate, the Downward/Northward/Eastward of the observation station
        // on or above the Earth's Surface,
        // Chnage the Polarity for: Gr, Gcolat, Grlon, Gcolatlon, Glonr, Gloncolat

        (*Gr)           = (-1.0) *K*( -ro + rs*costmp )/Ros3;
}

/*
    Given one observation point outside the model, compute the Gravity Effects from all the cells of model.
    Using the adaptive scheme according the distance from the source to the observation point.
*/
void OneGridGraGrAdaptive(MOD3DSphGra *mod3dsphgra, double ro, double lonrado, double colatrado, float *Gr )
{
    double  gr;
    double  rs1, lonrads1, colatrads1, rs2, lonrads2, colatrads2;
    double  rs, lonrads, colatrads;
    double  dlonrad, dcolatrad, dr;
    double  alonrad, acolatrad, ar;
    double   density;
    double   weight, wgl;

    int     nx, ny, nz;
    int     i, j, k;
    int     ii, jj, kk;

    /*
            parameters for Gauss-Legendre nodes number
    */
    double   mindis;
    float   dx, dy ,dz;
    int     xn, yn, zn;
    double  x1 = -1.0, x2 = 1.0;

    double  xw[NMAX+1], yw[NMAX+1], zw[NMAX+1];   //NOTICE:       xw, yw, zw, xp, yp, zp begin from 1 to xn, yn, zn
    double  xp[NMAX+1], yp[NMAX+1], zp[NMAX+1];   //NOTICE:       NOT from 0.
    k = NMAX+1 ;
    for(i=0; i<k; i++){
        xw[i] = 0.0; yw[i] = 0.0; zw[i] = 0.0;
        xp[i] = 0.0; yp[i] = 0.0; zp[i] = 0.0;
    }

    nx = mod3dsphgra->nx;
    ny = mod3dsphgra->ny;
    nz = mod3dsphgra->nz;

    (*Gr) = 0.0;

    for( k=(nz-1); k>=0; k-- ){//===from bottom to top==
    for( j=0; j<ny; j++ ){
    for( i=0; i<nx; i++ ){
        //get the beginning and ending points for each direction.
        if( k == (nz-1) ){
            rs1             = mod3dsphgra->dep[k];
            rs2             =0.5*(mod3dsphgra->dep[k-1]+mod3dsphgra->dep[k]);
        }
        else if ( k == 0 ){
            rs1             = 0.5*(mod3dsphgra->dep[k]+mod3dsphgra->dep[k+1]);
            rs2             = mod3dsphgra->dep[k];
        }
        else{
            rs1             = 0.5*(mod3dsphgra->dep[k]+mod3dsphgra->dep[k+1]);
            rs2             = 0.5*(mod3dsphgra->dep[k-1]+mod3dsphgra->dep[k]);
        }
                                                                                                                                                                                                
        if( i == 0 ){
            lonrads1        = mod3dsphgra->lon[i];
            lonrads2        = 0.5*(mod3dsphgra->lon[i+1]+mod3dsphgra->lon[i]);
        }
        else if( i == (nx-1) ){
            lonrads1        = 0.5* (mod3dsphgra->lon[i]+mod3dsphgra->lon[i-1]);
            lonrads2        = mod3dsphgra->lon[i];
        }
        else{
            lonrads1        = 0.5*(mod3dsphgra->lon[i]+mod3dsphgra->lon[i-1]);
            lonrads2        = 0.5*(mod3dsphgra->lon[i+1]+mod3dsphgra->lon[i]);
        }
                                                                                                                                                                                                
        if( j == 0 ){
            colatrads1      = mod3dsphgra->lat[j];
            colatrads2      = 0.5*(mod3dsphgra->lat[j+1]+mod3dsphgra->lat[j]);
        }
        else if ( j == (ny-1) ){
            colatrads1      = 0.5*(mod3dsphgra->lat[j]+mod3dsphgra->lat[j-1]);
            colatrads2      = mod3dsphgra->lat[j];
        }
        else{
            colatrads1      = 0.5*(mod3dsphgra->lat[j]+mod3dsphgra->lat[j-1]);
            colatrads2      = 0.5*(mod3dsphgra->lat[j+1]+mod3dsphgra->lat[j]);
        }
                                                                                                                                                                                                
        getdensityBlock( mod3dsphgra, i, j, k, &density );      //===get the density at the grid i/j/k block

        //====change to International Unit kg/m3
        //weight        = 1000.0*((rs2-rs1)*(lonrads2-lonrads1)*(colatrads2-colatrads1))/8.0;
        weight          = 125.0*((rs2-rs1)*(lonrads2-lonrads1)*(colatrads2-colatrads1));

        dlonrad         = lonrads2 - lonrads1;
        dcolatrad       = colatrads2 - colatrads1;
        dr              = rs2 - rs1;
        alonrad         = lonrads2 + lonrads1;
        acolatrad       = colatrads2 + colatrads1;
        ar              = rs2 + rs1;

        //====calculate the distance from source to observation point,
        //====determine the parameters for Gauss-Legendre quadrature integration
        dx = fabs( dlonrad*6371.0*sin(colatrads1) );
        dy = fabs( dcolatrad *6371.0 );
        dz = fabs( dr );

        cal_minDis( &mindis, colatrado, lonrado, ro, colatrads1, lonrads1, colatrads2, lonrads2, rs2 );
        if( mindis < MINDISKM )  mindis = MINDISKM;
        if(mindis<110.0 && density>1.0e-5){

            xn = (int) ( NSCALE/(mindis/dx) );//====set the scale of each small block is 1/10 of the distance to the observation point
            yn = (int) ( NSCALE/(mindis/dy) );
            zn = (int) ( NSCALE/(mindis/dz) );
            if( xn < NMIN ) xn = NMIN;
            if( yn < NMIN ) yn = NMIN;
            if( zn < NMIN ) zn = NMIN;
            if( xn > NMAX ) xn = NMAX;
            if( yn > NMAX ) yn = NMAX;
            if( zn > NMAX ) zn = NMAX;

            gauleg( x1, x2, xp, xw,  xn);
            gauleg( x1, x2, yp, yw,  yn);
            gauleg( x1, x2, zp, zw,  zn);

            //===do Gauss-Legendre quadrature integration===
            for( kk=1; kk<=zn; kk++ ){
            for( jj=1; jj<=yn; jj++ ){
            for( ii=1; ii<=xn; ii++ ){
                lonrads         = (xp[ii]*dlonrad + alonrad)/2.0;
                colatrads       = (yp[jj]*dcolatrad + acolatrad)/2.0;
                rs              = (zp[kk]*dr + ar )/2.0;
                wgl             = weight*xw[ii]*yw[jj]*zw[kk];

                //getdensity( mod3dsphgra, lonrads, colatrads, rs, &density );
                GravityR( ro, lonrado, colatrado, rs, lonrads, colatrads, density, &gr );
                (*Gr)  += wgl*gr;
            }}}
        }
    }}}

}

/*
    Given Random Observation System. Only Compute gravity in Radial direction.
*/
void OBSGravityGrRandomAdaptive(MOD3DSphGra *mod3dsphgra,OBSSphGraRandom *ObsSphGra )
{
    double   ro, lonrado, colatrado;
    int     i, j;
    int     np;
    //double  Gpotential, Gr, Gcolat, Glon, Grr, Grcolat, Grlon, Gcolatr, Gcolatcolat, Gcolatlon, Glonr, Gloncolat, Glonlon;

    if( ObsSphGra->israd==0 )
        chancoorOBSsphGraRandom( ObsSphGra, 1 );
    np  = ObsSphGra->np;

    printf("\nBegin calculating the gravity effects:\n");

    ro = ObsSphGra->z0;
    for( i=0; i<np; i++ ){

        lonrado   = ObsSphGra->lon[i];
        colatrado = ObsSphGra->lat[i];
        OneGridGraGrAdaptive( mod3dsphgra, ro, lonrado, colatrado, &ObsSphGra->Gr[i] );
        printf("computing %d-th point: lon=%f lat=%f,gr=%g\n",
              i+1,lonrado/degrad2,90-colatrado/degrad2,ObsSphGra->Gr[i]);
    }
}

/*
    Given one observation point outside the model, compute the Gravity matrix from all the cells of model.
    Using the adaptive scheme according the distance from the source to the observation point.
*/
void gravmat_one_row(MOD3DSphGra *mod3dsphgra, double ro, double lonrado,
                     double colatrado,int *col,float *rw,int *nar)
{
    double  gr;
    double  rs1, lonrads1, colatrads1, rs2, lonrads2, colatrads2;
    double  rs, lonrads, colatrads;
    double  dlonrad, dcolatrad, dr;
    double  alonrad, acolatrad, ar;
    double   density;
    double   weight, wgl,Gr;

    int     nx, ny, nz;
    int     i, j, k;
    int     ii, jj, kk;
    int n;

    /*
            parameters for Gauss-Legendre nodes number
    */
    double   mindis;
    float   dx, dy ,dz;
    int     xn, yn, zn;
    double  x1 = -1.0, x2 = 1.0;

    double  xw[NMAX+1], yw[NMAX+1], zw[NMAX+1];   //NOTICE:       xw, yw, zw, xp, yp, zp begin from 1 to xn, yn, zn
    double  xp[NMAX+1], yp[NMAX+1], zp[NMAX+1];   //NOTICE:       NOT from 0.
    k = NMAX+1 ;
    for(i=0; i<k; i++){
        xw[i] = 0.0; yw[i] = 0.0; zw[i] = 0.0;
        xp[i] = 0.0; yp[i] = 0.0; zp[i] = 0.0;
    }

    nx = mod3dsphgra->nx;
    ny = mod3dsphgra->ny;
    nz = mod3dsphgra->nz;

    for( k=(nz-1); k>=0; k-- ){//===from bottom to top==
    for( i=0; i<nx; i++ ){
    for( j=0; j<ny; j++ ){
        Gr = 0.0;
        n = k*nx*ny+i*ny+j;
        //get the beginning and ending points for each direction.
        if( k == (nz-1) ){
            rs1             = mod3dsphgra->dep[k];
            rs2             =0.5*(mod3dsphgra->dep[k-1]+mod3dsphgra->dep[k]);
        }
        else if ( k == 0 ){
            rs1             = 0.5*(mod3dsphgra->dep[k]+mod3dsphgra->dep[k+1]);
            rs2             = mod3dsphgra->dep[k];
        }
        else{
            rs1             = 0.5*(mod3dsphgra->dep[k]+mod3dsphgra->dep[k+1]);
            rs2             = 0.5*(mod3dsphgra->dep[k-1]+mod3dsphgra->dep[k]);
        }
                                                                                                                                                                                                
        if( i == 0 ){
            lonrads1        = mod3dsphgra->lon[i];
            lonrads2        = 0.5*(mod3dsphgra->lon[i+1]+mod3dsphgra->lon[i]);
        }
        else if( i == (nx-1) ){
            lonrads1        = 0.5* (mod3dsphgra->lon[i]+mod3dsphgra->lon[i-1]);
            lonrads2        = mod3dsphgra->lon[i];
        }
        else{
            lonrads1        = 0.5*(mod3dsphgra->lon[i]+mod3dsphgra->lon[i-1]);
            lonrads2        = 0.5*(mod3dsphgra->lon[i+1]+mod3dsphgra->lon[i]);
        }
                                                                                                                                                                                                
        if( j == 0 ){
            colatrads1      = mod3dsphgra->lat[j];
            colatrads2      = 0.5*(mod3dsphgra->lat[j+1]+mod3dsphgra->lat[j]);
        }
        else if ( j == (ny-1) ){
            colatrads1      = 0.5*(mod3dsphgra->lat[j]+mod3dsphgra->lat[j-1]);
            colatrads2      = mod3dsphgra->lat[j];
        }
        else{
            colatrads1      = 0.5*(mod3dsphgra->lat[j]+mod3dsphgra->lat[j-1]);
            colatrads2      = 0.5*(mod3dsphgra->lat[j+1]+mod3dsphgra->lat[j]);
        }
                                                                                                                                                                                                
        //getdensityBlock( mod3dsphgra, i, j, k, &density );      //===get the density at the grid i/j/k block

        //====change to International Unit kg/m3
        //weight        = 1000.0*((rs2-rs1)*(lonrads2-lonrads1)*(colatrads2-colatrads1))/8.0;
        weight          = 125.0*((rs2-rs1)*(lonrads2-lonrads1)*(colatrads2-colatrads1));

        dlonrad         = lonrads2 - lonrads1;
        dcolatrad       = colatrads2 - colatrads1;
        dr              = rs2 - rs1;
        alonrad         = lonrads2 + lonrads1;
        acolatrad       = colatrads2 + colatrads1;
        ar              = rs2 + rs1;

        //====calculate the distance from source to observation point,
        //====determine the parameters for Gauss-Legendre quadrature integration
        dx = fabs( dlonrad*6371.0*sin(colatrads1) );
        dy = fabs( dcolatrad *6371.0 );
        dz = fabs( dr );

        cal_minDis( &mindis, colatrado, lonrado, ro, colatrads1, lonrads1, colatrads2, lonrads2, rs2 );
        if( mindis < MINDISKM )  mindis = MINDISKM;
        if(mindis<100.0){

            xn = (int) ( NSCALE/(mindis/dx) );//====set the scale of each small block is 1/10 of the distance to the observation point
            yn = (int) ( NSCALE/(mindis/dy) );
            zn = (int) ( NSCALE/(mindis/dz) );
            if( xn < NMIN ) xn = NMIN;
            if( yn < NMIN ) yn = NMIN;
            if( zn < NMIN ) zn = NMIN;
            if( xn > NMAX ) xn = NMAX;
            if( yn > NMAX ) yn = NMAX;
            if( zn > NMAX ) zn = NMAX;

            gauleg( x1, x2, xp, xw,  xn);
            gauleg( x1, x2, yp, yw,  yn);
            gauleg( x1, x2, zp, zw,  zn);

            //===do Gauss-Legendre quadrature integration===
            for( kk=1; kk<=zn; kk++ ){
            for( jj=1; jj<=yn; jj++ ){
            for( ii=1; ii<=xn; ii++ ){
                lonrads         = (xp[ii]*dlonrad + alonrad)/2.0;
                colatrads       = (yp[jj]*dcolatrad + acolatrad)/2.0;
                rs              = (zp[kk]*dr + ar )/2.0;
                wgl             = weight*xw[ii]*yw[jj]*zw[kk];

                //getdensity( mod3dsphgra, lonrads, colatrads, rs, &density );
                GravityR( ro, lonrado, colatrado, rs, lonrads, colatrads, 1.0, &gr );
                Gr += wgl*gr;
            }}}
            if(Gr > 1.0e-6){
                rw[*nar] = Gr;
                col[*nar] = n;
                *nar = *nar + 1;
            }
        }
    }}}

}

/*
    Given Random Observation System. Only Compute matrix of gravity in Radial direction.
*/
int gravmat(MOD3DSphGra *mod3dsphgra,OBSSphGraRandom *ObsSphGra,float *rw,int *iw)
{
    double   ro, lonrado, colatrado;
    int     i, j;
    int     np,nar,nar1;
    int     *col;

    //double  Gpotential, Gr, Gcolat, Glon, Grr, Grcolat, Grlon, Gcolatr, Gcolatcolat, Gcolatlon, Glonr, Gloncolat, Glonlon;

    np  = ObsSphGra->np;
    ro = ObsSphGra->z0;
    nar = 0;
    nar1 = 0;
    col = ivec(np * (mod3dsphgra->nx * mod3dsphgra->ny * mod3dsphgra->nz));

    for( i = 0; i < np; i++ ){

        lonrado   = ObsSphGra->lon[i];
        colatrado = ObsSphGra->lat[i];
        gravmat_one_row(mod3dsphgra,ro,lonrado,colatrado,col,rw,&nar1);
        for(j = nar;j < nar1;j++){
            iw[j+1] = i;
        }
        printf("computing %d-th point: lon=%f lat=%f\n",
              i+1,lonrado/degrad2,90-colatrado/degrad2);
        nar = nar1;
    }
    iw[0] = nar;

    for(i = 0;i < nar;i++){
        iw[nar+1+i] = col[i];
    }

    free_ivec(col);

    return nar;
}

/*
    Given Random Observation System. Only Compute matrix of gravity in Radial direction.
    openmpi is used.
*/
int gravmat_parallel(MOD3DSphGra *mod3dsphgra,OBSSphGraRandom *ObsSphGra,float *rw,int *iw)
{
    double   ro;
    int     i, j;
    int     np,nar,nar1,m;

    np  = ObsSphGra->np;
    ro = ObsSphGra->z0;
    m = mod3dsphgra->nx * mod3dsphgra->ny * mod3dsphgra->nz;

    // allocate sparse matrix for every data
    int **col,*nars;
    float **data;

    col = imat(np,m);
    data = fmat(np,m);
    nars = ivec(np);

    #pragma omp parallel for
    for(i=0;i<np;i++){
        double lonrado   = ObsSphGra->lon[i];
        double colatrado = ObsSphGra->lat[i];
        gravmat_one_row(mod3dsphgra,ro,lonrado,colatrado,col[i],data[i],nars+i);
        printf("computing %d-th point: lon=%f lat=%f\n",
              i+1,lonrado/degrad2,90-colatrado/degrad2);
    }
    nar = 0,nar1 = 0;
    for(i=0;i<np;i++){
        nar1 += nars[i];
    }

    for(i=0;i<np;i++){
        for(j=0;j<nars[i];j++){
            iw[nar+1] = i;
            iw[nar+1+nar1] = col[i][j];
            rw[nar] = data[i][j];
            nar += 1;
        }
    }

    iw[0] = nar;

    free_imat(col,np);
    free_ivec(nars);
    free_fmat(data,np);
    return nar;

}

#undef	NMIN
#undef NMAX
#undef NSCALE
#undef MINDISKM
