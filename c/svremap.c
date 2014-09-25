#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "gaepsi.h"

int svremap_init(SVRemap * r, int remap[3][3]) {
    double BAS[3][3] = {0};
    int i, j;
    double uv, uu, uvp;
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            r->remap[i][j] = remap[i][j];
        } 
    }
    /* u1 = v1 */
    for(j = 0; j < 3; j++) {
        BAS[0][j] = remap[0][j];
    } 

    /* u2 = v2 - proj(u1, v2) */
    /* u3' = v3 - proj(u1, v3) */
    uv = 0;
    uvp = 0;
    uu = 0;
    for(j = 0; j < 3; j ++) {
        uv += remap[1][j] * BAS[0][j];
        uu += BAS[0][j] * BAS[0][j];
        uvp += remap[2][j] * BAS[0][j];
    }
    for(j = 0; j < 3; j++) {
        BAS[1][j] = remap[1][j] - uv / uu * BAS[0][j];
        BAS[2][j] = remap[2][j] - uvp / uu * BAS[0][j];
    } 
    /* u3 = u3' - proj(u2, u3p) */
    uv = 0;
    uu = 0;
    for(j = 0; j < 3; j ++) {
        uv += BAS[1][j] * BAS[2][j];
        uu += BAS[1][j] * BAS[1][j];
    }
    for(j = 0; j < 3; j++) {
        BAS[2][j] = BAS[2][j] - uv / uu * BAS[1][j];
    } 
    memset(r->T, 0, sizeof(double) * 16);
    for(i = 0; i < 3; i++) {
        uu = 0;
        for(j = 0; j < 3; j++) {
            uu += BAS[i][j] * BAS[i][j];
        }
        uu = sqrt(uu);
        for(j = 0; j < 3; j++) {
            r->T[i][j] = BAS[i][j] / uu;
        }
    }
    r->T[3][3] = 1.0;
    for(i = 0; i < 3; i++) {
        r->size[i] = 0;
        for(j = 0; j < 3; j++) {
            r->size[i] += remap[i][j] * r->T[i][j];
        }
        if (r->size[i] < 0) {
            /* make sure size is > 0 */
            r->size[i] *= -1;
            for(j = 0; j < 3; j++) {
                r->T[i][j] *= -1;
            }
        }
    }
    return 0;
}

/*
 * solve 0 < (I + x) * Q < b
 * x is [0, 1]
 *
 * I is the initial guess
 * y is new cooridnate
 * y = (I + x) * Q
 *
 * */

static double svremap_test(SVRemap * r, double x[3], double y[3], int I[3]) {
    double v1[4];
    double v2[4];
    int i;
    for(i = 0; i < 3; i ++) {
        v1[i] = x[i] + I[i];
    }
    v1[3] = 1.0;
    gtmat_apply(r->T, v1, v2);
    for(i = 0; i < 3; i ++) {
        y[i] = v2[i];
    }
    double badness = 0.0;
    for(i = 0; i < 3; i ++) {
        if(v2[i] < 0.0) {
            badness = fmax(badness, - v2[i]);
        }
        if(v2[i] > r->size[i]) {
            badness = fmax(badness, v2[i] - r->size[i]);
        }
    }
    return badness;
}
static void svremap_bounds(SVRemap * r, int Imin[3], int Imax[3]) {
    GTMatrix inv;
    gtmat_inverse(inv, r->T);
    int i;
    int a;
    for(i = 0; i < 27; i ++) {
        double v1[4];
        double v2[3];
        v1[3] = 1.0;
        for(a = 0; a < 3; a ++) {
            v1[a] = ((i >> a) & 1) * r->size[a];
        }
        gtmat_apply(inv, v1, v2);
        for(a = 0; a < 3; a ++) {
            if(floor(v2[a]) < Imin[a])  Imin[a] = floor(v2[a]);
            if(ceil(v2[a]) > Imax[a])  Imax[a] = ceil(v2[a]);
        }
    }
}
double svremap_apply(SVRemap * r, double x[3], double y[3], int I[3]) {
    if(svremap_test(r, x, y, I) == 0.0) return 0.0;
    int Imax[3] = {0}, Imin[3] = {0}, Ibest[3]; 
    double ybest[3];
    double bestbadness = 99.99; 
    int a;
    svremap_bounds(r, Imin, Imax);

    for(a = 0; a < 3; a++) {
        I[a] = Imin[a];
    }
    int done = 0;
    while(!done) {
        double badness = svremap_test(r, x, y, I);
        if(badness == 0.0) return 0.0;
        if(badness < bestbadness) {
            for(a = 0; a < 3; a++) {
                Ibest[a] = I[a];
                ybest[a] = y[a];
            }
            bestbadness = badness;
        }

        I[0] ++;
        for(a = 0; a < 3; a++) {
            if(I[a] != Imax[a] + 1) continue;
            if(a + 1 == 3) done = 1;
            I[a + 1] ++;
            I[a] = Imin[a];
        }
    }
    for(a = 0; a < 3; a++) {
        I[a] = Ibest[a];
        y[a] = ybest[a];
    }
    return bestbadness;
}    
#if 0
int main(int argc, char * argv[]) {
    int rm[3][3] = {
        13, -13, 12, 
        11, 13, 13, 
        4, 7, 5};
    /*
    int rm[3][3] = {
        1, 1, 0, 
        0, 1, 0, 
        0, 0, 1};*/
    SVRemap r;
    SVRemap_init(&r, rm);
    GTMatrix inv;
    GTMatrix check;
    gtmat_inverse(inv, r.T);
    Tmul(check, inv, r.T);
    //printf("%s\n", Tprint(r.T));
    //printf("check unitary %s\n", Tprint(check));
    printf("#new size %g %g %g\n", r.size[0], r.size[1], r.size[2]);
    double x, y, z;
    int I[3] = {0};
    for(x = 0; x <= 1.0; x += 0.1) {
    for(y = 0; y <= 1.0; y += 0.1) {
    for(z = 0; z <= 1.0; z += 0.1) {
        double v[3] = {x, y, z};
        double u[3];
        double badness = svremap_apply(&r, v, u, I);
        printf("%g %g %g %g %g %g %g\n", x, y, z, u[0], u[1], u[2], badness);
    }
    }
    }
    return 0;
}
#endif
