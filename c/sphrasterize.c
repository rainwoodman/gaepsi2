#include <stdio.h>
#include <math.h>
#include <string.h>
#include "gaepsi.h"

inline void gsph_image_write_single(GSPHImage * img, int x, int y, double * value) {
    int c;
    char * ptr = img->data;
    ptr += x * img->strides[0] + y * img->strides[1];
           
    for(c = 0; c < img->size[2]; c++, ptr += img->strides[2]) {
#pragma omp atomic
        *((float*) ptr) += value[c];
    }
}

inline void gsph_image_write_double(GSPHImage * img, int x, int y, double * value) {
    int c;
    char * ptr = img->data;
    ptr += x * img->strides[0] + y * img->strides[1];
           
    for(c = 0; c < img->size[2]; c++, ptr += img->strides[2]) {
#pragma omp atomic
        *((double*) ptr) += value[c];
    }
}

inline void gsph_image_write(GSPHImage * img, int x, int y, double * value) 
{
    if (img->itemsize == 8) {
        gsph_image_write_double(img, x, y, value);
    } else {
        gsph_image_write_single(img, x, y, value);
    }
}

void 
gsph_image_init(GSPHImage * image, 
        char * dtype, 
        int size[3], 
        ptrdiff_t * strides, 
        void * data) 
{
    if(dtype[strlen(dtype) - 1] == '8') {
        image->itemsize = 8;
    } else {
        image->itemsize = 4;
    }
    image->size[0] = size[0];
    image->size[1] = size[1];
    image->size[2] = size[2];
    if(strides) {
        image->strides[0] = strides[0];
        image->strides[1] = strides[1];
        image->strides[2] = strides[2];
    } else {
        /* assuming continuous strides */
        image->strides[2] = image->itemsize;
        image->strides[1] = image->strides[2] * size[2];
        image->strides[0] = image->strides[1] * size[1];
    }
    image->data = data;
}

void 
gsph_rasterize(GSPHImage * image, GSPHKernel sphkernel,
        double pos[2], double sml, double * mvalue) 
{
    /* sml here is half of Gadget's cubic spline sml (gadget uses support). */
    int * size = image->size;
    int nc = image->size[2];

    if (size[0] == 0) return;
    if (size[1] == 0) return;

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
        if (x >= size[1]) return;
        if (y >= size[0]) return;
        gsph_image_write(image, x, y, mvalue);
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
            if (min[k] >= size[k]) return;
            if (min[k] < 0) min[k] = 0;
            if (max[k] < 0) max[k] = 0;
            if (min[k] >= size[k]) min[k] = size[k] - 1;
            if (max[k] >= size[k]) max[k] = size[k] - 1;
        }
        if(sml < 60) {
            for(y = pos[0] - sml; y <= pos[0] + sml; y++) {
            for(x = pos[1] - sml; x <= pos[1] + sml; x++) {
                double dx = x - pos[1];
                double dy = y - pos[0];
                double r = sqrt(dx * dx + dy * dy) / (sml );
                r = sphkernel(r);
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
                double tmp[nc];
                int i;
                for(i = 0; i < nc; i ++) {
                    tmp[i] = mvalue[i] * w;
                }
                gsph_image_write(image, y, x, tmp);
            }
        }
        }
    }
}
