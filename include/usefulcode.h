double geocen2geogralat(double hlat);
double geogra2geocenlat(double hlat);
double compute_length(double colatErad, double lonErad, double RE, 
                    double colatSrad, double lonSrad, double RS );
void cal_minDis(double *mindis, double colatrado, double lonrado, double ro, 
            double colatrads1, double lonrads1, double colatrads2, double lonrads2, 
            double rs );

/*
    calloc array,matrix and 3-d tensor allocate function
*/
double *dvec(int n);
float *fvec(int n);
int *ivec(int n);
double **dmat(int row,int col);
float **fmat(int row,int col);
int **imat(int row,int col);
double ***d3tensor(int row,int col,int deep);

/*
    Deallocate function
*/
void free_dvec(double *a);
void free_fvec(float *a);
void free_ivec(int *a);
void free_dmat(double **P,int row_p);
void free_fmat(float **a,int row);
void free_imat(int **a,int row);
void free_d3tensor(double ***P,int row,int col);
