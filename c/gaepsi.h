#ifndef __GAEPSI_H__
#define __GAEPSI_H__
#include <stdint.h>
#include <stddef.h>
#define addressof(p) &(p)
/* geomtry and matrix operations, in matrix.c */
typedef double GTMatrix [4][4];

/* define a shift matrix */
#define GTMAT_SHIFT(x, y, z) { 1, 0, 0, x, 0, 1, 0, y, 0, 0, 1, z, 0, 0, 0, 1, }
/* define a boost matrix */
#define GTMAT_BOOST(x, y, z) { x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1, }

void gtmat_apply(GTMatrix T, double v1[4], double v2[4]);
void gtmat_mul(GTMatrix T, GTMatrix A, GTMatrix B);
void gtmat_copy(GTMatrix T, GTMatrix A);
void gtmat_inverse(GTMatrix T, GTMatrix A);

/* survey volume remapping, in svremap.c */
typedef struct {
    GTMatrix T;
    int remap[3][3];
    double size[3];
} SVRemap;

int svremap_init(SVRemap * r, int remap[3][3]);
double svremap_apply(SVRemap * r, double x[3], double y[3], int I[3]);


/* spline.c */
typedef double (*GSPHKernel)(double rr);
GSPHKernel * gsph_kernel_query(char * name);

/* sphrasterize.c */
typedef struct {
    /* multi channel image of float or double */
    int itemsize;
    /* x, y, c */
    int size[3]; 
    ptrdiff_t strides[3];
    void * data;
} GSPHImage;

void
gsph_image_init(GSPHImage * image, 
        char * dtype, 
        int size[3], 
        ptrdiff_t * strides, 
        void * data);

void 
gsph_rasterize(GSPHImage * image, GSPHKernel sphkernel,
        double pos[2], double sml, 
        double * mvalue) ;

/* quadcurve.c */
typedef int64_t quadindex_t;

quadindex_t quadindex(double x, double y);
#endif
