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
#define SPLINE_2D_PROJ_CUBIC "2dcubicproj"
typedef double (*gsph_spline_kernel)(double r);
gsph_spline_kernel gsph_spline_query(char * name);

/* sphrasterize.c */
typedef struct {
    int size[2]; 
    ptrdiff_t strides[2];
    char * wdata;
    char * vdata;
    int usedouble;
} GSPHImage;

typedef void (*gsph_painter_write)(void * data, int x, int y, double * values);

typedef struct {
    void * data;
    int size[2];
    gsph_painter_write write;
    gsph_spline_kernel sphkernel;
    int nvalue;
} GSPHPainter;

void gsph_image_init(GSPHImage * image, int size[2], char * dtype, 
        ptrdiff_t * strides, void * wdata, void * vdata);
void gsph_image_write(GSPHImage * image, int x, int y, double * mvalue);

void gsph_painter_init(GSPHPainter * painter, int size[2], 
        gsph_painter_write write, 
        gsph_spline_kernel sphkernel, 
        void * data, 
        int nvalue);

void gsph_painter_rasterize(GSPHPainter * painter,
        double pos[2], double sml, double * mvalue);

/* quadcurve.c */
typedef int64_t quadindex_t;

quadindex_t quadindex(double x, double y);
#endif
