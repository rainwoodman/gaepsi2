#include <stdio.h>
#include <math.h>
#include <string.h>
#include "gaepsi.h"

void gsph_painter_init(GSPHPainter * painter, int size[2], 
        gsph_painter_write write, 
        gsph_spline_kernel sphkernel, 
        void * data,
        int nvalue) {
    painter->size[0] = size[0];
    painter->size[1] = size[1];
    painter->write = write;
    painter->data = data;
    painter->sphkernel = sphkernel;
    painter->nvalue = nvalue;
}

void gsph_image_init(GSPHImage * image, int size[2], char * dtype, 
        ptrdiff_t * strides, void * wdata, void * vdata) {
    if(dtype[strlen(dtype) - 1] == '8') {
        image->usedouble = 1;
    } else {
        image->usedouble = 0;
    }
    image->size[0] = size[0];
    image->size[1] = size[1];
    if(strides) {
        image->strides[0] = strides[0];
        image->strides[1] = strides[1];
    } else {
        /* assuming continous strides */
        image->strides[1] = image->usedouble?sizeof(double):sizeof(float);
        image->strides[0] = image->strides[1] * size[1];
    }
    image->wdata = wdata;
    image->vdata = vdata;
}

void gsph_image_write(GSPHImage * image, int x, int y, double * mvalue) {
    ptrdiff_t p = y * image->strides[0] + x * image->strides[1];
    float * fpw = (float*) &image->wdata[p];
    float * fpv = (float*) &image->vdata[p];
#pragma omp atomic
    *fpw += mvalue[0];
#pragma omp atomic
    *fpv += mvalue[1];
}

void gsph_painter_rasterize(GSPHPainter * painter, 
        double pos[2], double sml, double * mvalue) {

    if (painter->size[0] == 0) return;
    if (painter->size[1] == 0) return;

    double bit = 0;
    int x, y;
    double r;
    int k;

    if(sml < 1.0) {
        /* really shall use CIC here fix it later */
        /* if sml is too small (less than 0.707) bits becomes zero */
        x = pos[1];
        y = pos[0];
        if (x < 0) return;
        if (y < 0) return;
        if (x >= painter->size[1]) return;
        if (y >= painter->size[0]) return;
        painter->write(painter->data, x, y, mvalue);
    } else {
        double save[128 * 128];
        int s = 0;
        int usekernel = 0;

        int min[2]; 
        int max[2]; 
        for(k = 0; k < 2; k ++) {
            min[k] = pos[k] - sml;
            max[k] = pos[k] + sml;
            if (max[k] < 0) return;
            if (min[k] >= painter->size[k]) return;
            if (min[k] < 0) min[k] = 0;
            if (max[k] < 0) max[k] = 0;
            if (min[k] >= painter->size[k]) min[k] = painter->size[k] - 1;
            if (max[k] >= painter->size[k]) max[k] = painter->size[k] - 1;
        }

        if(sml < 32) {
            for(y = pos[0] - sml; y <= pos[0] + sml; y++) {
            for(x = pos[1] - sml; x <= pos[1] + sml; x++) {
                double dx = x - pos[1];
                double dy = y - pos[0];
                r = sqrt(dx * dx + dy * dy) / sml;
                r = painter->sphkernel(r);
                bit += r;
                if(x >= min[1] && x <= max[1] 
                        && y >= min[0] && y <= max[0]) {
                    save[s] = r;
                    s++;
                }
            }
            }
            usekernel = 1;
        } else {
            /* particle to big, use a rectangular kernel */
            bit = (4 * sml * sml);
        }

        bit = 1.0 / bit;
        double sml2 = sml * sml;
        //printf("area %g %g %g %g\n", sml, bit, sml2, bit / sml2);
        s = 0;
        for(y = min[0]; y <= max[0]; y++) {
        for(x = min[1]; x <= max[1]; x++) {
            double w;
            if(usekernel) {
                w = save[s] * bit;
            } else {
                w = 1.0 * bit;
            }
            s ++;
            if(w) {
                double tmp[painter->nvalue];
                int i;
                for(i = 0; i < painter->nvalue; i ++) {
                    tmp[i] = mvalue[i] * w;
                }
                painter->write(painter->data, x, y, tmp);
            }
        }
        }
    }
}
