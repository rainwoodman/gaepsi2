#include <math.h>
#include <string.h>

#include "gaepsi.h"
#include <stdio.h>

static double gaussian_kernel(double rr);

GSPHKernel gsph_spline_query(char * name) {
    return gaussian_kernel;
}

/* The projected 2D cubic spline */
static double 
_gaussian_kernel(double r) 
{
    return exp(- 0.5 * (r * r * 4) );
}
static double _cubic_spline_2D_proj(double r) {
    double h2 = r * r;
    double h4 = h2 * h2;
    double h;
    double h3;
    const double fac = 1.0 / 1.00023;
    if (h2 < 0.25)  {
        //   1.909859 = 6 / pi
        return fac * (1.909859317102744 -10.23669021 * h2 -23.27182034 * h4 *(h2 - 1));
    }
    if (h2 < 1.0)  {
        h = sqrt(h2);
        h3 = h2 * h;
        //# -1/9. * 7/ pi * (1-x) **4 + 100 / 7. / pi * (1-x) **3 - 1.5 / pi * (1-x) ** 2 + 6 / 17. / 5 / pi* (1-x)
        return fac * (-0.2475743559207261 * h4 - 3.556986664656963 * h3 
                + 11.67894130021956 * h2 
                - 11.71909411592771* h 
                + 3.84471383628584);
    }
    return 0.0;
}
static double gaussian_kernel(double r) {
    static double stable[1025];
    static double * table = NULL;
    if(table == NULL) {
        int i;
        for(i = 0; i <= 1024; i ++) {
            double r = 1.0 * i / 1024;    
//            stable[i] = _gaussian_kernel(r);
            stable[i] = _cubic_spline_2D_proj(r);
        }
        table = stable;
    }
    if(r < 0) return 0;
    if(r >= 1) return 0;
    double i = r * 1024;
    int left = i;
    int right = left + 1;
    double u = i - left;
    double v = 1 - u;
    return table[left] * v + table[right] * u;
}
